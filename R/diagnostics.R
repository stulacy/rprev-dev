#' Estimate consistency of incidence data with an homogeneous Poisson process.
#'
#' This function compares the actual variance of the yearly incidence rates with rates
#' simulated from a Poisson process with overall rate equal to the overall mean rate.
#'
#' @param data Vector of absolute incidence for each complete year of the registry.
#' @param N_sim Integer to specify number of simulations.
#' @return Vector of p values for over- and under-dispersion based on the position of the
#' observed sequence variance in the distribution.
#' @examples
#' sim_check(incidence(registry_data$entrydate))
sim_check <- function(data, N_sim = 100000){
  var_sim <- vapply(seq(N_sim), function(i) var(rpois(length(data), mean(data))), numeric(1))
  c(length(var_sim[var_sim > var(data)])/N_sim, length(var_sim[var_sim <= var(data)])/N_sim)
}

#' Plot the age distribution of incident cases in the registry data.
#'
#' @param agedata Vector of age at diagnosis for each patient in the registry.
#' @param df Degrees of freedom for the smooth.
#' @return Plot of the raw data and smoothed age distribution function.
#' @examples
#' incidence_age_distribution(registry_data_r$age)
incidence_age_distribution <- function(agedata, df=10){

  ages <- vapply(seq(100), function(i) sum(floor(agedata) + 1 == i), numeric(1))

  plot(0:99, ages[1:100], pch=20, xlab="age (years)", ylab="incident cases")
  smage <- smooth.spline(0:99, ages[1:100], df=df)
  lines(smage, col="blue", lwd=2)

}

#' Inspect functional form of age.
#'
#' @param form ...
#' @param data A registry dataset of patient cases.
#' @param df Degrees of freedom for the smooth.
#' @return Plots of the functional form of age.
#' @examples
#' functional_form_age(registry_data_r)
functional_form_age <- function(form, data, df=4){

  ### TO DO/discuss:
  # ?No reason why this can't be applied to any continuous covariate, just need to change age() and age_ prefixes
  # ?How to neaten up the output; control side effects, do we need both plots etc
  # ?Too much duplication of code here with prevalence()
  ###

  # Extract required column names from formula
  terms <- terms(form, 'age')
  special_indices <- attr(terms, 'specials')

  if (any(sapply(special_indices, is.null)))
    stop("Error: provide function term for age.")

  v <- as.list(attr(terms, 'variables'))[-1]
  var_names <- unlist(lapply(special_indices, function(i) v[i]))
  age_var <- .extract_var_name(var_names$age)

  # Extract survival formula
  response_index <- attr(terms, 'response')
  resp <- v[response_index][[1]]

  psp_surv_form <- as.formula(paste(deparse(resp), '~ pspline(',
                                    age_var, ', ', df, ')', sep=''))

  cxnl <- coxph(psp_surv_form, data)
  output1 <- summary(cxnl)

  plt1 <- termplot(cxnl)

  f <<- datadist(data)
  options(datadist="f")

  rcs_surv_form <- as.formula(paste(deparse(resp), '~ rcs(',
                                    age_var, ', ', df, ')', sep=''))

  mod_rms <- cph(rcs_surv_form, data, x=TRUE, y=TRUE, surv=T, time.inc=1)
  output2 <- anova(mod_rms)
  output3 <- summary(mod_rms)

  plt2 <- plot(eval(parse(text=paste('Predict(mod_rms, ', age_var,')', sep = ''))), lwd=3, adj.subtitle=T)

  list(plt1, plt2, output1, output2, output3)

}

#' ?
#'
#' @param age Age of example patient.
#' @param sex Sex of example patient ("Male" or "Female").
#' @param wb Weibull model fitted to the data.
#' @return ?
dfr_sc <- function(age, sex, wb){
  if(sex == "Male"){
    sexn <- 0
  }else{
    sexn <- 1
  }
  return(exp(wb$coef[1] + age*wb$coef[2] + sexn*wb$coef[3]))
}

#' ?
#'
#' @param t ?
#' @param shape ?
#' @param scale ?
#' @return ?
wb_Su <- function(t, shape, scale){
  return(exp(-(t/scale)^shape))
}

#' ?
#'
#' @param t ?
#' @param shape ?
#' @param scale ?
#' @return ?
wb_haz <- function(t, shape, scale){
  return((shape/scale)*(t/scale)^(shape-1))
}

#' ?
#'
#' @param age Age of example patient.
#' @param sex Sex of example patient ("Male" or "Female").
#' @param the_length Number of days to predict.
#' @param wb Weibull model fitted to the data.
#' @return ?
plot_haz <- function(age, sex, the_length=3000, wb){
  s_time <- seq(0, the_length, by=1)
  haz <- wb_haz(s_time, dfr_sh, dfr_sc(age, sex, wb))
  plot(s_time, haz, type="l", lwd=2, col="red", ylim=c(0,0.002),
       main=paste("age = ", age, ", sex =", sex), xlab="survival (days)", ylab="probability")
}

#' ?
#'
#' @param data Population survival dataset.
#' @param age Age of example patient.
#' @param sex Sex of example patient ("Male" or "Female").
#' @param the_length Number of days to predict.
#' @return ?
plot_population_hazard <- function(data, age, sex, the_length=3000){
  if(sex == "Male") sex = "Males"
  if(sex == "Female") sex = "Females"
  daily_survival_rate <- function(data, sex)
  lines(0:(the_length-1), daily_survival_rate[(365*age):(365*age + the_length - 1)],
        col="blue", lwd=2)
}

#' ?
#'
#' @param age Age of example patient.
#' @param sex Sex of example patient ("Male" or "Female").
#' @param the_length Number of days to predict.
#' @param wb Weibull model fitted to the data.
#' @return ?
plot_Su <- function(age, sex, the_length=3000, wb){
  s_time <- seq(0, the_length, by=1)
  Su <- wb_Su(s_time, dfr_sh, dfr_sc(age, sex, wb))
  plot(s_time, Su, type="l", lwd=2, col="red", ylim=c(0,1),
       main=paste("age = ", age, ", sex =", sex), xlab="survival (days)", ylab="probability")
}

#' ?
#'
#' @param data Population survival dataset.
#' @param age Age of example patient.
#' @param sex Sex of example patient ("Male" or "Female").
#' @param the_length Number of days to predict.
#' @return ?
plot_pop_Su <- function(data, age, sex, the_length=3000){
  if(sex == "Male") sex = "Males"
  if(sex == "Female") sex = "Females"
  daily_survival_rate <- function(data, sex)
  lines(1:the_length, cumprod(1 - daily_survival_rate[(age*365):(age*365 + the_length - 1)]),
        col="cyan", lwd=2)
}

#' ?
#'
#' @param data A registry dataset of patient cases generated using load_data().
#' @param registry_years A vector of dates delineating years of the registry.
#' @param registry_start_year Ordinal defining the first year of the registry data to be used.
#' @param age Age of example patient.
#' @param sex Sex of example patient ("Male" or "Female").
#' @param limits ?
#' @param colour ?
#' @return ?
plot_km <- function(data, registry_years, registry_start_year, age,
                    sex, limits, colour="green"){
  if(sex == "Male"){
    sexn <- 0
  }else{
    sexn <- 1
  }
  dfr_r <- data[data$date_initial >= registry_years[registry_start_year],]
  kma <- survfit(Surv(survival_time, indicator) ~ 1,
                 data=dfr_r[dfr_r$sex == sexn & dfr_r$age_initial > age-limits &
                              dfr_r$age_initial < age + limits, ])
  lines(kma, lwd=1, col=colour, conf.int=T)
}

#' ?
#'
#' @param data A registry dataset of patient cases generated using load_data().
#' @param registry_years A vector of dates delineating years of the registry.
#' @param registry_start_year Ordinal defining the first year of the registry data to be used.
#' @param age Age of example patient.
#' @param sex Sex of example patient ("Male" or "Female").
#' @param N_boot Number of replicates.
#' @return ?
boot_eg <- function(data, registry_years, registry_start_year, age, sex, N_boot = 1000){

  wb_boot <- registry_survival_bootstrapped(data = data[data$date_initial >= registry_years[registry_start_year], ])
  wb_lines <- matrix(0, nrow=N_boot, ncol=5000)
  n <- seq(1,5000, by=1)

  if(sex == "Male"){
    for (i in 1:N_boot){
      wb_lines[i, ] <- 1 - pweibull(n, scale=exp(wb_boot[i, 1] + age*wb_boot[i, 2]), shape=1/wb_boot[i, 4])
    }
  } else {
    for (i in 1:N_boot){
      wb_lines[i, ] <- 1 - pweibull(n, scale=exp(wb_boot[i, 1] + age*wb_boot[i, 2] + wb_boot[i, 3]), shape=1/wb_boot[i, 4])
    }
  }

  plot(NA, xlim=c(1,5000), ylim=c(0,1), main = paste("age = ", age, ", sex = ", sex))
  for (i in 1:N_boot){
    lines(wb_lines[i, ], lwd=1, col="grey")
  }

  wb <- survreg(Surv(survival_time, indicator) ~ age_initial + sex, data)

  if(sex == "Male"){
    A <- 1 - pweibull(n, scale=exp(wb$coef[1] + age*wb$coef[2]), shape=1/wb$scale)
    lines(A, col="orange", lwd=3, lty=3)
  } else {
    A <- 1 - pweibull(n, scale=exp(wb$coef[1] + age*wb$coef[2] + wb$coef[3]), shape=1/wb$scale)
    lines(A, col="orange", lwd=3, lty=3)
  }

}

#' Chi squared test between prevalence prediction and observed values in the registry.
#'
#' @param entry Vector of diagnosis dates for each patient in the registry.
#' @param events Vector of event or censorship dates for each patient in the registry.
#' @param status Vector of binary values indicating if an event has occurred for each patient in the registry.
#' @param by_year Vector of predicted number of prevalent cases by each year of diagnosis.
#' @param start Date from which incident cases are included.
#' @param num_reg_years Integer representing the number of complete years of the registry for which incidence is to be calculated.
#' @return Chi-squared test of difference between prevalence prediction and counted prevalence at the index date.
#' @examples
#' prev_chisq(entry = registry_data$entrydate,
#'            events = registry_data$eventdate,
#'            status = registry_data$status,
#'            start="2004-01-30", num_reg_years = 9,
#'            by_year = by_year_total)
prev_chisq <- function(entry, events, status, by_year, start=NULL, num_reg_years=NULL){
  observed <- counted_prevalence(entry, events, status, start, num_reg_years)
  predicted <- rev(by_year[1:num_reg_years])
  chi <- sum(((observed - predicted)^2)/predicted)
  1 - pchisq(chi, num_reg_years - 1)

}