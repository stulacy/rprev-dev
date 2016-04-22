#' Count prevalence from registry data.
#'
#' @param data A registry dataset of patient cases generated using load_data().
#' @param registry_years A vector of dates delineating years of the registry.
#' @param registry_start_year Ordinal defining the first year of the registry data to be used.
#' @param registry_end_year Ordinal  defining the last year of the registry data to be used.
#' @return A count of prevalence at the index date subdivided by year of diagnosis and inclusion in the registry.
#' @examples
#' counted_prevalence(load_data(registry_data), registry_years, registry_start_year, registry_end_year)
counted_prevalence<-function(data, registry_years, registry_start_year, registry_end_year){

  data <- data[!(is.na(data$date_event)), ]
  events <- sum(data$indicator)
  if (events < 30) warning("warning: low number of events.")

  years_estimated <- registry_end_year - registry_start_year + 1
  per_year <- rep(NA, years_estimated)

  for (i in registry_start_year:registry_end_year){
    per_year[i - registry_start_year + 1] <-
      length(data$date_initial[data$date_initial >= registry_years[i] & data$date_initial < registry_years[i + 1]])
  }

  counted_prevalence <- rep(NA, years_estimated)
  registry_years <- registry_years[registry_start_year:(registry_start_year + years_estimated)]

  for (i in 1:years_estimated){
    counted_prevalence[i] <- per_year[i] -
      sum(data$indicator_censored_at_index[data$date_initial >= registry_years[i] & data$date_initial < registry_years[i + 1]])
  }

  return(counted_prevalence)

}

#' Predict prevalence for a given number of years.
#'
#' @param data A registry dataset of patient cases generated using load_data().
#' @param registry_years A vector of dates delineating years of the registry.
#' @param registry_start_year Ordinal defining the first year of the registry data to be used.
#' @param registry_end_year Ordinal defining the last year of the registry data to be used.
#' @param daily_survival_males Population survival function (Males).
#' @param daily_survival_females Population survival function (Females).
#' @param cure_time Cure model assumption for the calculation.
#' @param age_division Can be "None", "Adult" or "Paediatric".
#' @param age_cut A number, in years, to specify the cut-off between adult and
#'        paediatric cases.
#' @param sex Can be "Both", "Males" or "Females".
#' @param N_years Number of years to model.
#' @param N_boot Number of bootstrapped calculations to perform.
#' @param Max_Yearly_Incidence A number larger than the expected yearly incidence
#'        to allow for variation in incidence between years.
#' @return Predicted number of prevalent cases for each year, the age distribution of the
#'        prevalent cases, an indication of sampling variation and raw incidence data to
#'        allow for crosschecking with previous calculations.
#' @examples
#' by_year_total <- prevalence(load_data(registry_data), registry_years, registry_start_year,
#'                            registry_end_year, daily_survival_males = daily_survival_males,
#'                            daily_survival_females = daily_survival_females, cure_time = cure*365,
#'                            N_years = 10)
prevalence <- function(data, registry_years, registry_start_year, registry_end_year, N_years,
                       daily_survival_males, daily_survival_females, cure_time,
                       N_boot = 1000, Max_Yearly_Incidence = 500){

  if(mean(data$sex) == 0 | mean(data$sex) == 1) type = "single"
  if(mean(data$sex) != 0 | mean(data$sex) != 1) type = "both"

  data_r <- data[data$date_initial >= registry_years[registry_start_year], ]
  est_years <- registry_end_year - registry_start_year + 1

  if(cure_time > 0){
    wb_boot <- registry_survival_bootstrapped(data = data[data$date_initial >= registry_years[registry_start_year], ], N_boot)
  }else{
    wb <- survreg(Surv(survival_time, indicator) ~ age_initial + sex, data=data)
    wb_boot <- t(vapply(1:N_boot, function(x) c(wb$coe, wb$scale), FUN.VALUE = numeric(4)))
  }
  wb_boot <- wb_boot[sample(nrow(wb_boot)), ]

  if(type != "both"){

    fix_rate <- incidence(data, registry_years, registry_start_year, registry_end_year)
    if(mean(data$sex) == 0){
      prior_age_d <- data_r$age_initial[data_r$sex == 0]
    }else{
      prior_age_d <- data_r$age_initial[data_r$sex == 1]
    }
    post_age_dist <- array(NA, dim = c(N_years, N_boot, Max_Yearly_Incidence))
    by_year_samples <- matrix(NA, nrow=N_years, ncol=N_boot)
    fix_rate_rev <- rev(fix_rate)
    mean_rate <- mean(fix_rate)

    for (j in 1:N_years){
      year_no <- j - 1
      boot_out <- rep(NA, N_boot)
      for (i in 1:N_boot){
        rate <- max(0, rnorm(1, mean_rate, sqrt(mean_rate)/est_years))
        if(j <= est_years){
          rate <- fix_rate_rev[j]
        }

        no_diag <- rpois(1, rate)
        boot_age_dist <- sample(prior_age_d, no_diag, replace=T)

        diag_time <- year_no * 365 + runif(no_diag, 0, 365)

        if(mean(data$sex) == 0){
          d_or_a <- rbinom(no_diag, 1, 1 - survival_model(diag_time, boot_age_dist, sex = 0,
                                                   cure_time, boot = wb_boot[i,], daily_survival = daily_survival_males))
        }else{
          d_or_a <- rbinom(no_diag, 1, 1 - survival_model(diag_time, boot_age_dist, sex = 1,
                                                     cure_time, boot = wb_boot[i,], daily_survival = daily_survival_females))
        }

        K_ages <- length(diag_time[d_or_a == 0])
        if(K_ages > 0) post_age_dist[j , i, 1:K_ages] <-
          diag_time[d_or_a == 0]/365 + boot_age_dist[d_or_a == 0]
        boot_out[i] <- no_diag - sum(d_or_a)
      }

      by_year_samples[j, ] <- boot_out

    }

  }else{

    prior_age_d_male <- data_r$age_initial[data_r$sex == 0]
    prior_age_d_female <- data_r$age_initial[data_r$sex == 1]

    post_age_dist_male <- array(NA, dim = c(N_years, N_boot, Max_Yearly_Incidence))
    post_age_dist_female <- array(NA, dim = c(N_years, N_boot, Max_Yearly_Incidence))

    by_year_male_samples <- matrix(NA, nrow=N_years, ncol=N_boot)
    by_year_female_samples <- matrix(NA, nrow=N_years, ncol=N_boot)

    fix_m_rate <- rev(incidence(data = data[data$sex == 0, ], registry_years, registry_start_year, registry_end_year))
    fix_f_rate <- rev(incidence(data = data[data$sex == 1, ], registry_years, registry_start_year, registry_end_year))

    mean_rate_male <- mean(fix_m_rate)
    mean_rate_female <- mean(fix_f_rate)

    for (j in 1:N_years){
      year_no <- j - 1
      boot_out_male <- rep(NA, N_boot)
      boot_out_female <- rep(NA, N_boot)
      for (i in 1:N_boot){
        male_rate <- max(0, rnorm(1, mean_rate_male, sqrt(mean_rate_male)/est_years))
        female_rate <- max(0, rnorm(1, mean_rate_female, sqrt(mean_rate_female)/est_years))
        if(j <= (registry_end_year - registry_start_year + 1)){
          male_rate <- fix_m_rate[j]
          female_rate <- fix_f_rate[j]
        }
        no_diag_male <- rpois(1, male_rate)
        no_diag_female <- rpois(1, female_rate)
        boot_age_dist_male <- sample(prior_age_d_male, no_diag_male, replace=T)
        boot_age_dist_female <- sample(prior_age_d_female, no_diag_female, replace=T)

        diag_time_male <- year_no * 365 + runif(no_diag_male, 0, 365)
        diag_time_female <- year_no * 365 + runif(no_diag_female, 0, 365)
        d_or_a_male <- rbinom(no_diag_male, 1,
                              1 - survival_model(diag_time_male, boot_age_dist_male, sex = 0,
                                          cure_time, boot = wb_boot[i,], daily_survival = daily_survival_males))
        d_or_a_female <- rbinom(no_diag_female, 1,
                                1 - survival_model(diag_time_female, boot_age_dist_female, sex = 1,
                                              cure_time, boot = wb_boot[i,], daily_survival = daily_survival_females))
        K_ages_m <- length(diag_time_male[d_or_a_male == 0])
        K_ages_f <- length(diag_time_female[d_or_a_female == 0])
        if(K_ages_m > 0) post_age_dist_male[j , i, 1:K_ages_m] <-
          diag_time_male[d_or_a_male == 0]/365 + boot_age_dist_male[d_or_a_male == 0]
        if(K_ages_f > 0) post_age_dist_female[j , i, 1:K_ages_f] <-
          diag_time_female[d_or_a_female == 0]/365 + boot_age_dist_female[d_or_a_female == 0]
        boot_out_male[i] <- no_diag_male - sum(d_or_a_male)
        boot_out_female[i] <- no_diag_female - sum(d_or_a_female)
      }
      by_year_male_samples[j, ] <- boot_out_male
      by_year_female_samples[j, ] <- boot_out_female

    }

    by_year_samples <- by_year_male_samples + by_year_female_samples
    post_age_dist <- abind(post_age_dist_male, post_age_dist_female, along=2)
    fix_rate <- rev(fix_m_rate + fix_f_rate)

  }

  by_year <- apply(by_year_samples, 1, mean)

  return(list(by_year, post_age_dist, by_year_samples, fix_rate))

}

#' Calculate predicted prevalence for a given number of years.
#'
#' @param N_years Number of years prevalence is to be calculated for.
#' @param registry_start_year Ordinal defining the first year of the registry data to be used.
#' @param registry_end_year Ordinal defining the last year of the registry data to be used.
#' @param the_samples Bootstrapped samples used to gauge the sampling variation in the estimates.
#' @param by_year Predicted prevalence by year of the registry.
#' @param population_size A number corresponding to the size of the population.
#' @return An estimate of prevalence, in total and /100,000 with confidence intervals.
#' @examples
#' n_year_estimates(N_years = 3, registry_start_year = year1, registry_end_year = lastyear,
#'                  the_samples = by_year_male_samples,
#'                  by_year = by_year_male, population_size = 1700000)
n_year_estimates <- function(N_years, registry_start_year, registry_end_year,
                             the_samples, by_year, population_size){
  if (N_years > length(by_year)) stop("error: too many years for the data.")

  the_samples <- the_samples[(registry_end_year - registry_start_year + 2):N_years, , drop=F]
  by_sample_estimate <- apply(the_samples, 2, sum)

  the_estimate <- sum(by_year[1:N_years])

  raw_proportion <- the_estimate / population_size
  the_proportion <- 100000 * raw_proportion

  if (N_years <= (registry_end_year - registry_start_year + 1)){
    CI <- (1.96 * sqrt((raw_proportion * (1 - raw_proportion))/population_size)) * 100000
  }else {
    the_estimate_n <- sum(by_year[1:(registry_end_year - registry_start_year + 1)])
    raw_proportion_n <- the_estimate_n / population_size
    std_err_1 <- sqrt((raw_proportion_n * (1 - raw_proportion_n))/population_size)
    std_err_2 <- sd(by_sample_estimate)/population_size

    CI <<- 1.96 * sqrt(std_err_1^2 + std_err_2^2) * 100000
  }

  raw <- paste("raw: ", round(the_estimate, 2), sep="")
  prop <- paste("/100,000 cases: ", round(the_proportion, 2),
        " (", round(the_proportion - CI, 2),
        "-", round(the_proportion + CI, 2), ")", sep="")

  return(list(raw, prop))
}

#' Output extrapolated number of prevalent cases per 10-year age group.
#'
#' @param dist The age distribution of prevalent cases following simulation.
#' @param registry_end_year The last complete year of registry data.
#' @param N_years The year to make predictions until.
#' @return Predicted number of prevalent cases per 10-year age group.
#' @examples
#' result <- prevalence_by_age(dist = post_age_dist_male, registry_end_year = 4,
#'                  N_years = N_years)
prevalence_by_age <- function(dist, registry_end_year, N_years){

  the_dist <- dist[registry_end_year:N_years, , ]
  the_dist <- unlist(the_dist)
  the_dist <- the_dist[!is.na(the_dist)]

  aL <- length(the_dist)

  a_0_9 <- length(the_dist[the_dist > 0 & the_dist < 10])/aL
  a_10_19 <- length(the_dist[the_dist >= 10 & the_dist < 20])/aL
  a_20_29 <- length(the_dist[the_dist >= 20 & the_dist < 30])/aL
  a_30_39 <- length(the_dist[the_dist >= 30 & the_dist < 40])/aL
  a_40_49 <- length(the_dist[the_dist >= 40 & the_dist < 50])/aL
  a_50_59 <- length(the_dist[the_dist >= 50 & the_dist < 60])/aL
  a_60_69 <- length(the_dist[the_dist >= 60 & the_dist < 70])/aL
  a_70_79 <- length(the_dist[the_dist >= 70 & the_dist < 80])/aL
  a_80_plus <- length(the_dist[the_dist >= 80])/aL

  return(c(a_0_9,  a_10_19, a_20_29, a_30_39, a_40_49, a_50_59,
              a_60_69, a_70_79, a_80_plus))

}
