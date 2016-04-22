#' Synthetic patient dataset.
#'
#' A dataset synthesised to resemble a disease registry where patient cases are
#' represented by a unique identifier. Demographic and disease-specific data
#' required for prevalence estimations are included, such as sex, age, and dates
#' of diagnosis and death or censorship. Deceased patients have \code{status = 1}, and the
#' entry in \code{EventDate} is the date of their death. Patients censored at some point
#' have \code{status = 0}, and the entry in \code{EventDate} is the date of their
#' censorship. The counted prevalence calculation requires a specified index date after
#' which all patient cases who are not deceased are censored, in this case 31st August 2013.
#' The corresponding recoded indicator variable is \code{status2}.
#'
#' @format A data frame with 1908 rows and 8 variables:
#' \describe{
#'  \item{ID}{ID, a unique identifier for the patient}
#'  \item{sex}{sex, "0" for males or "1" for females}
#'  \item{age}{age, in years}
#'  \item{DateOfDiag}{date of diagnosis}
#'  \item{status}{binary variable, 1 if patient is deceased and 0 if alive or censored}
#'  \item{EventDate}{date of death or censorship}
#'  \item{stime}{time between date of diagnosis and death or censorship, in days}
#'  \item{status2}{numeric variable, 1 if deceased and 0 if censored at index date, otherwise NA}
#' }
"registry_data"

#' Load registry data.
#'
#' @param data A registry dataset of patient cases. Cases with missing values will be removed.
#' @param age_cut A number to define age cut-off in adult or paediatric analyses.
#' @param age_initial The column name representing age at data_initial.
#' @param date_initial The column name representing the date cases are added to the registry.
#' @param indicator The column name representing a binary variable where 1 is event has occurred and 0 is event has not occurred, or the case is censored.
#' @param survival_time The column representing time to event or censorship in days.
#' @param date_event The column name representing the date an event or censorship has occurred.
#' @param indicator_censored_at_index The column name representing a binary variable where 1 is event has occurred and 0 is event has not occurred, or the case is censored, by an index date of interest.
#' @param sex Can be "Both" for all cases, "Males" or "Females". Please note that the column name must be exactly "sex" to enable correct filtering.
#' @param age_division Can be "None", "Adult" or "Paediatric".
#' @return A dataset restricted to a population of desired age and sex. Selected age, sex, date and indicator columns are named appropriately for use in subsequent functions within this package.
#' @examples
#' load("../data/registry_data.rda")
#' summary(load_data(registry_data))
load_data <- function(data,
                      age_initial = "age",
                      date_initial = "DateOfDiag",
                      indicator = "status",
                      survival_time = "stime",
                      date_event = "EventDate",
                      indicator_censored_at_index = "status2",
                      sex = "Both",
                      age_division = "None",
                      age_cut = 20){

  if(!sex %in% c('Both', 'Males', 'Females')) stop("error: incorrect sex.")
  if(sex == "Males") data = filter(data, sex == 0)
  if( sex == "Females") data = filter(data, sex == 1)

  if(!age_division %in% c('None', 'Adult', 'Paediatric')) stop("error: incorrect age_division.")
  if(age_division == "Adult") data = filter(data, age >= age_cut)
  if(age_division == "Paediatric") data = filter(data, age < age_cut)

  myCols <- c(age_initial, "sex", date_initial, indicator, survival_time, date_event, indicator_censored_at_index)
  data <- data %>%
    select(match(myCols, names(.)))

  names(data) <- c("age_initial",
                   "sex",
                   "date_initial",
                   "indicator",
                   "survival_time",
                   "date_event",
                   "indicator_censored_at_index")
  return(data)
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
#' incidence(load_data(registry_data), registry_years, registry_start_year = 1, registry_end_year = 4)
incidence_dev <- function(entry, registry_years) {
  per_year <- vapply(seq(length(registry_years)-1),
                     function(i) sum(entry >= registry_years[i] & entry < registry_years[i+1]),
                     integer(1))

  if(sum(per_year) < 30) warning("warning: low number of incident cases.")
  return(per_year)
}

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
#' @param data A registry dataset of patient cases generated using load_data().
#' @param registry_years A vector of dates delineating years of the registry.
#' @param registry_start_year Ordinal defining the first year of the registry data to be used.
#' @param registry_end_year Ordinal defining the last year of the registry data to be used.
#' @param population_size A number representing the size of the population at risk.
#' @param precision The number of decimal places required.
#' @return A vector containing mean incidence/year and incidence rate/100,000/year with confidence
#'         intervals.
#' @examples
#' incidence_rates <- meanIR(load_data(registry_data), registry_years, registry_start_year = 1,
#' registry_end_year = 4, population_size = 3500000)
meanIR_dev <- function(entry, registry_years, population_size, precision = 2, level=0.95){

  mean_rate <- mean(incidence_dev(entry, registry_years))
  z_conf <- qnorm((1+level)/2)
  
  CI <- 100000 * (z_conf * sqrt(mean_rate) / length(registry_years)) / population_size
  est <- 100000 * mean_rate / population_size

  object <- list(absolute=mean_rate, per100K=est, per100K.lower=est-CI, per100K.upper=est+CI)
  lapply(object, round, precision)
}

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