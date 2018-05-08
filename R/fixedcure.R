#' Fixed cure model where survival reverts to population mortality
#' at a specified timepoint.
#'
#' @inheritParams flexsurvcure_population
#' @param formula Survival formula
#' @param data Dataset
#' @param dist Distribution to be used
#' @param time Time-limit at which a patient is considered cured. Uses
#'     same time-scale as survival time in \code{data}.
#'
#' @return An object of class \code{fixedcure} that can be passed
#'   into \code{\link{prevalence}}.
fixed_cure <- function(formula, data, dist, time, pop_data) {
    obj <- build_survreg(formula, data, user_dist)

    # TODO Make a function to validate mortality data so can reuse it with flexsurvcure
    if (is.null(pop_data)) {
        utils::data('UKmortalitydays', envir=environment())
        pop_data <- get('UKmortalitydays', envir=environment())
    }

    # Save individual attributes (NOT AGE) that are in mortality data
    # TODO This should identify all covariates that are in pop_data that AREN'T 'surv'
    mortality_covars <- setdiff(intersect(obj$covars, colnames(pop_data)), 'age')
    obj$pop_covars <- mortality_covars

    num_rates <- pop_data %>%
                    dplyr::count_(mortality_covars)
    if (any(num_rates$n != 36525))
        stop("Error: 'pop_data' must have 36525 rows for daily mortality over 100 years for each covariate combination")

    obj$pop_data <- pop_data
    obj$cure_time <- time

    class(obj) <- 'fixedcure'

    obj
}

# TODO Need consistent time scale, i.e. YEARS or DAYS! Ideally wouldn't matter, but does make a difference for the cure model
# Maybe just have the documentation state that cure time needs to be on same time scale as that used in Surv(time, status)

#' @export
predict_survival_probability.fixedcure <- function(object, newdata, times) {

    raw_probs <- predict_survival_probability.survregmin(object, newdata, times)

    # TODO Find times which are greater than cure time

    # Calculate population mortalities at these times

}

#' @export
extract_covars.fixedcure <- function(object) {
    object$covars
}