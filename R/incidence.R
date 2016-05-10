#' Synthetic patient dataset.
#'
#' A dataset synthesised to resemble a disease registry of patient cases. Demographic and
#' disease-specific variables required for prevalence estimation are included, such as sex, age,
#' and dates of diagnosis, death and the corresponding indicator.
#'
#' @format A data frame with 1000 rows and 6 variables:
#' \describe{
#'  \item{time}{time between date of diagnosis and death or censorship, in days}
#'  \item{status}{binary variable, 1 if patient is deceased and 0 if alive or censored}
#'  \item{age}{age, in years}
#'  \item{sex}{sex, "0" for males and "1" for females}
#'  \item{entrydate}{date of diagnosis}
#'  \item{eventdate}{date of death or censorship, corresponding to \code{status}}
#' }
"prevsim"

#' Calculate absolute incidence from registry data.
#'
#' @param entry Vector of diagnosis dates for each patient in the registry.
#' @param start Date from which incident cases are included.
#' @param num_reg_years Integer representing the number of complete years of the registry for which incidence is to be calculated.
#' @return Vector of absolute incidence values for each included year of the registry.
#' @examples
#' incidence(registry_data$DateOfDiag, start="2005-09-01", 8)
incidence <- function(entry, start=NULL, num_reg_years=NULL) {

    if (is.null(start))
        start <- min(entry)

    if (is.null(num_reg_years))
        num_reg_years <- floor(as.numeric(difftime(max(entry), start) / 365.25))

    registry_years <- .determine_registry_years(start, num_reg_years)

    per_year <- vapply(seq(num_reg_years),
                       function(i) sum(entry >= registry_years[i] & entry < registry_years[i+1]),
                       integer(1))

    if(sum(per_year) < 30) warning("warning: low number of incident cases.")
    per_year
}

#' Calculate absolute incidence from registry data.
#'
#' @param data A registry dataset of patient cases generated using load_data().
#' @param registry_years A vector of dates delineating years of registry.
#' @param registry_start_year Ordinal defining the first year of registry data to be used.
#' @param registry_end_year Ordinal defining the last year of registry data to be used.
#' @return Vector of absolute incidence values for each of the specified years of a registry.
#' @examples
#' registry_years <- c("2004-09-01", "2005-09-01", "2006-09-01", "2007-09-01", "2008-09-01")
#' incidence_current(load_data(registry_data), registry_years, registry_start_year = 1, registry_end_year = 4)
incidence_current <- function(data, registry_years, registry_start_year, registry_end_year){
    years_estimated <- registry_end_year - registry_start_year + 1
    per_year <- rep(NA, years_estimated)

    for (i in registry_start_year:registry_end_year){
        per_year[i - registry_start_year + 1] <-
            length(data$date_initial[data$date_initial >= registry_years[i] & data$date_initial < registry_years[i + 1]])
    }
    total <- sum(per_year)
    if(total < 30) warning("warning: low number of incident cases.")
    return(per_year)
}

#' Convert absolute incidence values to per 100,000 population values.
#'
#' @param entry Vector of diagnosis dates for each patient in the registry.
#' @param population_size Integer representing the size of the population at risk.
#' @param start Date from which incident cases are included.
#' @param num_reg_years Integer representing the number of complete years of the registry for which incidence is to be calculated.
#' @param precision Integer representing the number of decimal places required.
#' @param level Double representing the desired confidence interval width.
#' @return List containing mean incidence/year and incidence rate/100,000/year with confidence intervals.
#' @examples
#' incidence_rates <- meanIR(load_data(registry_data), registry_years, registry_start_year = 1,
#' registry_end_year = 4, population_size = 3500000)
meanIR <- function(entry, population_size, start=NULL, num_reg_years=NULL, precision = 2, level=0.95){

    if (is.null(start))
        start <- min(entry)

    if (is.null(num_reg_years))
        num_reg_years <- floor(as.numeric(difftime(max(entry), start) / 365.25))

    registry_years <- .determine_registry_years(start, num_reg_years)

    mean_rate <- mean(incidence(entry, start, num_reg_years=num_reg_years))
    z_conf <- qnorm((1+level)/2)

    CI <- 100000 * (z_conf * sqrt(mean_rate) / num_reg_years) / population_size
    est <- 100000 * mean_rate / population_size

    object <- list(absolute=mean_rate, per100K=est, per100K.lower=est-CI, per100K.upper=est+CI)
    lapply(object, round, precision)
}

#' Convert absolute incidence values to per 100,000 population values.
#'
#' @param data A registry dataset of patient cases generated using load_data().
#' @param registry_years A vector of dates delineating years of the registry.
#' @param registry_start_year Ordinal defining the first year of the registry data to be used.
#' @param registry_end_year Ordinal defining the last year of the registry data to be used.
#' @param population_size A number representing the size of the population at risk.
#' @param precision The number of decimal places required.
#' @return A vector containing mean incidence/year and incidence rate/100,000/year with confidence
#'         intervals.
#' @examples
#' incidence_rates <- meanIR_current(load_data(registry_data), registry_years, registry_start_year = 1,
#' registry_end_year = 4, population_size = 3500000)
meanIR_current <- function(data, registry_years, registry_start_year, registry_end_year, population_size, precision = 2){

    per_year <- incidence_current(data, registry_years, registry_start_year, registry_end_year)
    mean_rate <- mean(per_year)
    CI <- 100000 * (1.96 * sqrt(mean_rate) / (registry_end_year - registry_start_year + 1)) / population_size
    est <- 100000 * mean_rate / population_size

    result <- data.frame(round(c(mean_rate, est - CI, est, est + CI),
                               precision))
    rownames(result) <- c("absolute", "per100000-CI", "per100000", "per100000+CI")
    return(result)

}