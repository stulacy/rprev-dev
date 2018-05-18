DAYS_IN_YEAR <- 365.25

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

# Checks that all the factors to stratify population survival by are in the registry data set.
validate_population_survival <- function(population_data, registry_data, population_covariates=NULL) {

    if (!'surv' %in% colnames(population_data))
        stop("Error: no 'surv' column found in population survival rates 'daily_survival'.")

    if (min(population_data$surv) < 0 | max(population_data$surv) > 1)
        stop("Error: population survival probabilities found in 'surv' are out of the range [0,1].")

    if (!'age' %in% colnames(population_data))
        stop("Error: no 'age' column found in population survival rates 'daily_survival'.")
    if (!'age' %in% colnames(registry_data))
        stop("Error: no 'age' column found in registry data frame 'data'.")

    # Remove age from population covariates as must be there anyway
    population_covariates <- setdiff(population_covariates, 'age')

    pop_covariates_in_pop <- sapply(population_covariates, function(x) x %in% colnames(population_data))
    pop_covariates_in_reg <- sapply(population_covariates, function(x) x %in% colnames(registry_data))

    if (!all(pop_covariates_in_pop))
        stop(paste0("Error: not all values in population_covariates are in 'daily_survival'. Missing ",
                   paste(population_covariates[!pop_covariates_in_pop], collapse=","), "."))

    if (!all(pop_covariates_in_reg))
        stop(paste0("Error: not all values in population_covariates are in 'data'. Missing ",
                   paste(population_covariates[!pop_covariates_in_reg], collapse=","), "."))

    missing_covar_factors <- sapply(population_covariates, function(covar) {
        all(unique(registry_data[[covar]]) %in% population_data[[covar]])
    })
    
    if (!all(missing_covar_factors))
        stop(paste("Error: not all values of", 
                   paste(population_covariates[!missing_covar_factors], collapse=','),
                   "are present in population survival."))
}