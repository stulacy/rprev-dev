#' Builds a cure model with an associated population mortality table
#'
#' @param pop_data A data frame comprising population survival as a daily rate for 36525 days (100 years).
#'   It can also be stratified by other variables that are found in the survival \code{formula} for this model,
#'   such as sex.
#' @inheritParams flexsurvcure::flexsurvcure
#'
#' @importFrom magrittr "%>%"
#' @export
flexsurvcure_population <- function(formula, data, dist='weibull', pop_data=NULL, link='logistic', mixture=T, ...) {
    obj <- flexsurvcure::flexsurvcure(formula=formula, data=data,
                                      dist=dist,
                                      link=link, mixture=mixture,
                                      ...)

    if (is.null(pop_data)) {
        utils::data('UKmortalitydays', envir=environment())
        pop_data <- get('UKmortalitydays', envir=environment())
    }

    # Save individual attributes (NOT AGE) that are in mortality data

    model_covs <- attr(obj$concat.formula, "covnames")
    mortality_covars <- setdiff(intersect(model_covs, colnames(pop_data)), 'age')
    obj$pop_covars <- mortality_covars

    num_rates <- pop_data %>%
                    dplyr::count_(mortality_covars)
    if (any(num_rates$n != 36525))
        stop("Error: 'pop_data' must have 36525 rows for daily mortality over 100 years for each covariate combination")

    # Save pop_data to obj
    obj$pop_mortality <- pop_data
    obj$call <- match.call()

    obj
}

#' @export
predict_survival_probability.flexsurvcure <- function(object, newdata=NULL,
                                                      times=NULL)
{
    # Get estimates from flexsurvreg method
    estimates <- predict_survival_probability.flexsurvreg(object, newdata, times)

    if (!is.null(object$pop_mortality)) {
        # Obtain age at index in days
        # TODO Possible to not have hardcoded?
        newdata$age <- floor((newdata$age * 365.25) + t)

        # Obtain mortality rates from these values
        comb <- suppressWarnings(dplyr::left_join(newdata, object$pop_mortality, by=c('age', object$pop_covars)))

        # For people that don't have a mortality set to 0 as assumedly because they
        # are > 100
        comb[is.na(comb$surv), 'surv'] <- 0

        # scale estimates
        estimates <- estimates * comb$surv
    }

    estimates
}
