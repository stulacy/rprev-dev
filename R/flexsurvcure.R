#' @description Fits mixture and non-mixture cure models using \code{flexsurvcure}. Please read the detailed description below
#'   for how to use this model.
#' @inherit curemodels title details params
#' @param formula Formula specifying survival model as used by \code{\link{flexsurvcure}{flexsurvcure}}.
#' @inheritParams flexsurvcure::flexsurvcure
#' @return An object of class \code{flexsurvcure}.
#' @importFrom magrittr "%>%"
#' @export
flexsurvcure_population <- function(formula, data, daily_survival=NULL, population_covariates=NULL,
                                    dist='weibull', link='logistic', mixture=TRUE, ...) {
    obj <- flexsurvcure::flexsurvcure(formula=formula, data=data,
                                      dist=dist,
                                      link=link, mixture=mixture,
                                      ...)

    if (is.null(daily_survival)) {
        utils::data('UKmortalitydays', envir=environment())
        daily_survival <- get('UKmortalitydays', envir=environment())
    }

    validate_population_survival(daily_survival, data, population_covariates)

    obj$pop_covars <- population_covariates
    obj$pop_mortality <- daily_survival
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
        newdata$age <- floor((newdata$age * DAYS_IN_YEAR) + times)

        # Obtain mortality rates from these values
        comb <- suppressWarnings(dplyr::left_join(newdata, object$pop_mortality, by=c('age', object$pop_covars)))

        # For people that don't have a mortality set to 0 as assumedly because they
        # are > 100
        # TODO more explicit solution for this please!
        comb[is.na(comb$surv), 'surv'] <- 0

        # scale estimates
        estimates <- estimates * comb$surv
    }

    estimates
}
