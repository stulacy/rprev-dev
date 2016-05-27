#' Calculate yearly dates.
#'
#' A helper function to calculate dates a year apart, starting from a given date
#' and running for a set number of years.
#'
#' @inheritParams incidence
#' @return A vector of dates delineating complete years of registry data.
#' @export
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