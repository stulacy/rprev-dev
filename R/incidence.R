#' Simulated patient dataset.
#'
#' A dataset in the format of a disease registry, where the outcome being modeled is death
#' due to the simulated disease. The registry began in March 2013, with 1000 incident cases being
#' recorded over a period of nearly exactly ten years. The patients are followed up for a further two years
#' until 17.03.2015, at which point any subjects alive are marked as right censored.
#'
#' Demographic and disease-specific data required for prevalence estimations are included,
#' such as sex, age, and dates of entry and event. \code{eventdate} marks the date of the
#' last known follow-up with the patient, which is either death (\code{status = 1}) or
#' right-censored (\code{status = 0}).
#'
#' @format A data frame with 1000 rows and 6 columns:
#' \describe{
#'  \item{time}{(double) time between date of diagnosis and death or censorship, in days}
#'  \item{status}{(integer) event marker, 1 if patient is deceased and 0 if alive or censored}
#'  \item{age}{(double) age in years at point of entry into the registry}
#'  \item{sex}{(integer), 0 for males or 1 for females}
#'  \item{entrydate}{(character) date of entry into the registry in YYYY-MM-DD format}
#'  \item{eventdate}{(character) date of death or censorship in YYYY-MM-DD format}
#' }
"prevsim"

#' Calculate absolute incidence from registry data.
#'
#' @param entry (character): Vector of diagnosis dates for each patient in the registry in the format YYYY-MM-DD.
#' @param start (character): Date from which incident cases are included in the format YYYY-MM-DD, defaults to the earliest entry date.
#' @param num_reg_years (integer): The number of years of the registry for which incidence is to be calculated, defaults to using all available complete years.
#' @return Vector of length num_reg_years of doubles, representing the number of absolute incidence values for each included year of the registry.
#' @examples
#' data(prevsim)
#'
#' incidence(prevsim$entrydate, start="2004-01-01", 8)
#' incidence(prevsim$entrydate)
#' incidence(prevsim$entrydate, start="2005-05-01", 5)
#' incidence(prevsim$entrydate, start="2005-05-01")
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

#' Convert absolute incidence values to per 100,000 population values.
#'
#' @param entry (character): Vector of diagnosis dates for each patient in the registry in the format YYYY-MM-DD.
#' @param population_size (integer): The size of the population at risk.
#' @param start (character): Date from which incident cases are included in the format YYYY-MM-DD, defaults to the earliest entry date.
#' @param num_reg_years (integer): The number of years of the registry for which incidence is to be calculated, defaults to using all available complete years.
#' @param precision (integer): The number of decimal places required.
#' @param level (double): The desired confidence interval width.
#' @return A list with the following values:
#'
#' \item{absolute}{Overall incidence for the period of interest as a single double}
#' \item{per100K}{Incidence for the period of interest per one hundred thousand}
#' \item{per100K.lower}{Lower bounds of the specified confidence level on the per one hundred thousand estimate}
#' \item{per100K.upper}{Upper bounds of the specified confidence level on the per one hundred thousand estimate}
#'
#' @examples
#' data(prevsim)
#'
#' meanIR(prevsim$entrydate, population_size=3.5e6, start='2004-01-01', num_reg_years=8)
#' meanIR(prevsim$entrydate, population_size=3.5e6)
#' meanIR(prevsim$entrydate, population_size=3.5e6, precision=3)
#' meanIR(prevsim$entrydate, population_size=3.5e6, start='2004-01-01', num_reg_years=8, level=0.99)
meanIR <- function(entry, population_size, start=NULL, num_reg_years=NULL, precision = 2, level=0.95){

    if (is.null(start))
        start <- min(entry)

    if (is.null(num_reg_years))
        num_reg_years <- floor(as.numeric(difftime(max(entry), start) / 365.25))

    registry_years <- .determine_registry_years(start, num_reg_years)

    mean_rate <- mean(incidence(entry, start, num_reg_years=num_reg_years))
    z_conf <- qnorm((1+level)/2)

    CI <- 100e3 * (z_conf * sqrt(mean_rate) / num_reg_years) / population_size
    est <- 100e3 * mean_rate / population_size

    object <- list(absolute=mean_rate, per100K=est, per100K.lower=est-CI, per100K.upper=est+CI)
    lapply(object, round, precision)
}