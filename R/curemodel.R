#' Builds a cure model with an associated population mortality table
#'
#' @param pop_data A data frame comprising population survival as a daily rate for 36525 days (100 years).
#'   It can also be stratified by other variables that are found in the survival \code{formula} for this model,
#'   such as sex.
#' @inheritParams flexsurvcure::flexsurvcure
#'
#' @importFrom magrittr "%>%"
#' @export
cure_model_wrapper <- function(formula, data, dist='weibull', pop_data=NULL, link='logistic', mixture=T, ...) {
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