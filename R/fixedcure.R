#' @description Fits a cure model which assumes that if an individual has survived
#'  beyond a set time-point then they are considered cured and their mortality reverts
#'  to population levels.
#'  Please read the detailed description below for how to use this model.
#' @inherit curemodels title details params
#' @inheritParams prevalence
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
    func_call$dist <- eval(dist)
    func_call$time <- eval(time)
    func_call$daily_survival <- eval(daily_survival)
    obj$call <- func_call

    if (is.null(daily_survival)) {
        utils::data('UKmortalitydays', envir=environment())
        daily_survival <- get('UKmortalitydays', envir=environment())
    }

    validate_population_survival(daily_survival, data, population_covariates)

    # Save individual attributes (NOT AGE) that are in mortality data
    # TODO This should identify all covariates that are in daily_survival that AREN'T 'surv'
    mortality_covars <- setdiff(intersect(obj$covars, colnames(daily_survival)), 'age')
    obj$pop_covars <- mortality_covars

    obj$pop_data <- daily_survival
    obj$cure_time <- time

    class(obj) <- 'fixedcure'

    obj
}

#' @export
predict_survival_probability.fixedcure <- function(object, newdata, times) {

    browser()
    raw_probs <- predict_survival_probability.survregmin(object, newdata, times)
    pop_indices <- setNames(object$pop_covars, object$pop_covars)

    cured_individuals <- newdata %>%
        dplyr::mutate(time_to_index = times,
                      raw_prob = raw_probs) %>%
        dplyr::filter(time_to_index > object$cure_time)

    if (nrow(cured_individuals) > 0) {

        cured_individuals <- cured_individuals %>%
            dplyr::mutate(age_at_cure = floor(age * DAYS_IN_YEAR + object$cure_time),
                          age_at_index = floor(age * DAYS_IN_YEAR + time_to_index))
        # TODO How to get this working for people who have age at cure or age at index > max age in mortality data? Assume dead?
        adj_probs <- cured_individuals %>%
                dplyr::left_join(object$pop_data, by=c('age_at_cure'='age', pop_indices)) %>%
                dplyr::left_join(object$pop_data, by=c('age_at_index'='age', pop_indices), suffix=c('.cure', '.index')) %>%
                dplyr::mutate(adj_prob = raw_prob * (surv.index / surv.cure))

        raw_probs[times > object$cure_time] <- adj_probs$adj_prob
    }

    raw_probs
}

#' @export
extract_covars.fixedcure <- function(object) {
    object$covars
}