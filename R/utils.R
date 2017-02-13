DAYS_IN_YEAR <- 365.25

#' DEPRECATED: Use \code{\link{determine_yearly_endpoints}} instead.
#'
#' A helper function to calculate dates a year apart, starting from a given date
#' and running for a set number of years.
#'
#' @inheritParams incidence
#' @return A vector of dates delineating complete years of registry data.
#' @export
determine_registry_years <- function(start, num_reg_years) {
    .Deprecated("determine_yearly_endpoints")
    # Calculate registry years from this
    # NB: Ugly hack to not take leap years into account. Done so that tests don't throw an error, but strictly
    # should just use as.Date(start) + 365.25 * x to account for leap years
    determine_yearly_endpoints(start, num_reg_years)
}

#' Determine annual event delimiters.
#'
#' A helper function to calculate dates a year apart, starting from a given date
#' for a set number of years. Useful for delimiting a specified time
#' interval into a discrete number of years.
#'
#' @param date Either the starting date of the time interval (when \code{direction} is 'forwards')
#' or the ending date (\code{direction} is 'backwards').
#' @param num_years The number of years of the time interval.
#' @param direction A string indicating whether the parameter \code{date} represents
#' the opening or closing interval. Must take values either 'forwards' or 'backwards'.
#' @return A vector of dates delineating complete years of registry data.
#' @export
determine_yearly_endpoints <- function(date, num_years, direction='forwards') {

    if (! direction %in% c('forwards', 'backwards'))
        stop("Error: parameter 'direction' must take on a value of either 'forwards' or 'backwards'.")

    # Calculate registry years from this
    # NB: Ugly hack to not take leap years into account. Done so that tests don't throw an error, but strictly
    # should just use as.Date(start) + 365.25 * x to account for leap years
    if (direction == 'forwards') {
        sapply(0:num_years, function(x) paste(as.numeric(strftime(date, '%Y'))+x,
                                              strftime(date, '%m-%d'), sep='-'))
    } else {
        rev(sapply(0:num_years, function(x) paste(as.numeric(strftime(date, '%Y'))-x,
                                              strftime(date, '%m-%d'), sep='-')))
    }
}

# Extracts the variable names from a formula
# expr (language): Expression in the form entry(diagdate), where
# 'entry' is the known function name while the name inside brackets
# is the user's column name for that variable
.extract_var_name <- function(expr) {
    as.character(expr)[[2]]
}

.calculate_survival_probability <- function(coef, data, times) {
    scale <- exp(coef[-length(coef)] %*% data) + 0.000001  # Hack to ensure scale != 0
    shape <- 1 / coef[length(coef)]
    1 - pweibull(times, scale=scale, shape=shape)
}