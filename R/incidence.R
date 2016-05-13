#' Calculate absolute incidence from registry data.
#'
#' @param entry Vector of diagnosis dates for each patient in the registry in the format YYYY-MM-DD.
#' @param start Date from which incident cases are included in the format YYYY-MM-DD,
#' defaults to the earliest entry date.
#' @param num_reg_years The number of years of the registry for which incidence is to be calculated,
#' defaults to using all available complete years.
#' @return Vector of length num_reg_years of integers, representing the number of absolute incidence values
#' for each included year of the registry.
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

    registry_years <- determine_registry_years(start, num_reg_years)

    per_year <- vapply(seq(num_reg_years),
                       function(i) sum(entry >= registry_years[i] & entry < registry_years[i+1]),
                       integer(1))

    if(sum(per_year) < 30) warning("Warning: low number of incident cases.")
    per_year
}

#' Calculate average incidence rates per one hundred thousand over the specified number of years, with
#' confidence intervals.
#'
#' @param entry Vector of diagnosis dates for each patient in the registry in the format YYYY-MM-DD.
#' @param population_size The size of the population at risk.
#' @param start Date from which incident cases are included in the format YYYY-MM-DD, defaults to the earliest entry date.
#' @param num_reg_years The number of years of the registry for which incidence is to be calculated, defaults to using all available complete years.
#' @param precision The number of decimal places required.
#' @param level The desired confidence interval width.
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
#' mean_incidence_rate(prevsim$entrydate, population_size=3.5e6, start='2004-01-01', num_reg_years=8)
#' mean_incidence_rate(prevsim$entrydate, population_size=3.5e6)
#' mean_incidence_rate(prevsim$entrydate, population_size=3.5e6, precision=3)
#' mean_incidence_rate(prevsim$entrydate, population_size=3.5e6, start='2004-01-01', num_reg_years=8, level=0.99)
mean_incidence_rate <- function(entry, population_size, start=NULL, num_reg_years=NULL, precision = 2, level=0.95){

    if (is.null(start))
        start <- min(entry)

    if (is.null(num_reg_years))
        num_reg_years <- floor(as.numeric(difftime(max(entry), start) / 365.25))

    registry_years <- determine_registry_years(start, num_reg_years)

    mean_rate <- mean(incidence(entry, start, num_reg_years=num_reg_years))
    z_conf <- qnorm((1+level)/2)

    CI <- 100e3 * (z_conf * sqrt(mean_rate) / num_reg_years) / population_size
    est <- 100e3 * mean_rate / population_size

    object <- list(absolute=mean_rate, per100K=est, per100K.lower=est-CI, per100K.upper=est+CI)
    lapply(object, round, precision)
}