#' Builds survival models for diseases with cured fractions using population mortality tables
#'
#' Fits a cure model which assumes that if an individual has survived
#'  beyond a set time-point then they are considered cured and their mortality reverts
#'  to population levels.
#'  Please read the detailed description below for how to use this model.
#'
#' To model population survival, population mortality tables are required, as specified by the \code{daily_survival}
#' argument. If not provided, then the default population mortality is that of the UK population, which goes up
#' to 100 years of age. If a simulated individual has expected lifespan longer than the maximum age in the mortality table
#' then they are estimated to have died at this age limit,
#' which is why it is advantageous to provide as many accurate survival probabilities as possible.
#'
#' Due to the linking with the registry data and the ability for user-specified mortality tables, there are stricter
#' requirements on the survival models used in cure models than elsewhere. For example, the time-scale of the
#' survival model specified in \code{formula} \strong{must be in days} so that it matches up with the mortality tables.
#' Likewise, \strong{age in years must be included as a covariate} in the survival model
#'
#' @inheritParams prevalence
#' @param daily_survival A data frame comprising population survival as a daily probability for as long as possible,
#'   ideally 100 years (36525 days).
#'   Defaults to using UK mortality from the \code{UKmortalitydays} data set.
#'   It \strong{must contain columns 'age' and 'surv'}, providing the age (in days) and survival
#'   probability at that age respectively.
#'   It can also be stratified by other variables that are found in the survival \code{formula} for this model,
#'   such as sex.
#' @param population_covariates A character vector containing fields to stratify population survival by in addition to
#'   age, as descripted in \code{Details} below. These \strong{must} be the names of columns in both \code{data} and
#'   \code{daily_survival}. If not provided then defaults to the fields that are present in both \code{data} and
#'   \code{daily_survival}.
#' @param cure_time Time-limit at which a patient is considered cured.
#' @param formula Formula specifying survival function, as used in
#'   \code{\link{prevalence}} with the \code{surv_formula} argument.
#'   \strong{Must be in days}.
#'
#' @return An object of class \code{fixedcure} that can be passed
#'   into \code{\link{prevalence}}.
#' @export
fixed_cure <- function(formula, data, cure_time, daily_survival=NULL, population_covariates=NULL,
                       dist=c('exponential', 'weibull', 'lognormal')) {
    dist <- match.arg(dist)
    obj <- build_survreg(formula, data, dist)

    func_call <- match.call()
    func_call$formula <- eval(formula)
    obj$call <- func_call

    if (is.null(daily_survival)) {
        utils::data('UKmortalitydays', envir=environment())
        daily_survival <- get('UKmortalitydays', envir=environment())
    }

    # If don't provide population covariates then set to the ones in both data sets
    if (is.null(population_covariates)) {
        population_covariates <- setdiff(intersect(colnames(data), colnames(daily_survival)), 'age')
    }

    validate_population_survival(daily_survival, data, population_covariates)

    obj$pop_covars <- population_covariates
    obj$pop_data <- daily_survival
    obj$cure_time <- cure_time

    class(obj) <- 'fixedcure'

    obj
}

#' @export
predict_survival_probability.fixedcure <- function(object, newdata, times) {

    # For R CMD CHECK
    time_to_index <- NULL
    age <- NULL
    raw_prob <- NULL
    surv.index <- NULL
    surv.cure <- NULL
    surv_prob <- NULL
    adj_prob <- NULL
    time_calc_survival <- NULL
    id <- NULL
    cured <- NULL
    . <- NULL

    # Calculate raw survival time
    probs <- newdata %>%
        dplyr::mutate(id = 1:nrow(newdata),
                      time_to_index = times,
                      cured = time_to_index > object$cure_time,
                      time_calc_survival = ifelse(cured, object$cure_time, time_to_index),
                      surv_prob = predict_survival_probability.survregmin(object, ., time_calc_survival))

    if (sum(probs$cured) > 0) {
        # Cured survival probabilty is defined by Simon as
        # S(t) = S(curetime) * S*(t) / S*(curetime)
        # Where S*(.) is population survival rates calculated in terms of a person's age
        # S(curetime) was already calculated for these cured patients above and is held in 'surv_prob'

        pop_indices <- setNames(object$pop_covars, object$pop_covars)

        # Calculate the population survival rates by linking to the population survival rates
        cured_probs <- probs %>%
                    dplyr::filter(cured) %>%
                    dplyr::mutate(age_at_cure = floor(age * DAYS_IN_YEAR + object$cure_time),
                                  age_at_index = floor(age * DAYS_IN_YEAR + time_to_index))

        # Suppress warnings when joining by character field
        cured_probs_2 <- suppressWarnings(dplyr::left_join(cured_probs, object$pop_data, by=c('age_at_cure'='age', pop_indices)))
        cured_probs_3 <- suppressWarnings(dplyr::left_join(cured_probs_2, object$pop_data, by=c('age_at_index'='age', pop_indices),
                                                           suffix=c('.cure', '.index')))
        cured_probs_4 <- cured_probs_3 %>%
                            dplyr::mutate(adj_prob = surv_prob * (surv.index / surv.cure),
                                          adj_prob = ifelse(is.na(adj_prob), 0, adj_prob))

        # Link back with all newdata
        probs <- probs %>%
                    dplyr::left_join(cured_probs_4 %>% dplyr::select(id, adj_prob), by='id') %>%
                    dplyr::mutate(surv_prob = ifelse(cured, adj_prob, surv_prob))
    }

    # Return a vector
    probs %>% dplyr::pull(surv_prob)
}

#' @export
extract_covars.fixedcure <- function(object) {
    object$covars
}

print.fixedcure <- function(object, ...) {
    to_print <- object[c('coefs', 'covars', 'call', 'dist', 'terms', 'pop_covars', 'cure_time')]
    to_print$population_survival <- summary(object$pop_data)
    print(to_print)

}