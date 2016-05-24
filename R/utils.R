#' Find dates delineating complete years of registry data.
#'
#' @param start Date from which incident cases are included in the format YYYY-MM-DD.
#' Defaults to the earliest entry date.
#' @param num_reg_years The number of years of the registry for which incidence is to be calculated.
#' Defaults to using all available complete years.
#' @return A vector of dates delineating complete years of registry data.
#' @examples
#' determine_registry_years(start='2004-01-30', num_reg_years=8)
determine_registry_years <- function(start, num_reg_years) {
    # Calculate registry years from this
    # NB: Ugly hack to not take leap years into account. Done so that tests don't throw an error, but strictly
    # should just use as.Date(start) + 365.25 * x to account for leap years
    sapply(0:num_reg_years, function(x) paste(as.numeric(strftime(start, '%Y'))+x,
                                          strftime(start, '%m-%d'), sep='-'))
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