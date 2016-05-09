#' Estimate consistency of incidence data with a homogeneous Poisson process.
#'
#' This function compares the actual variance of the yearly incidence rates with rates simulated from a Poisson process with overall rate equal to the overall mean rate.
#'
#' @param data Vector of absolute incidence for each complete year of the registry.
#' @param N_sim Integer to specify number of simulations.
#' @return Vector of p values for over- and under-dispersion based on the position of the observed sequence variance in the distribution.
#' @examples
#' sim_check(raw_incidence)
sim_check <- function(data, N_sim = 100000){
  the_mean_rate <-  mean(data)
  var_sim <- rep(NA, N_sim)
  M <- length(data)
  for (i in 1:N_sim){
    thesims <- rpois(M, the_mean_rate)
    var_sim[i] <- var(thesims)
  }
  return(c(length(var_sim[var_sim > var(data)])/N_sim, length(var_sim[var_sim <= var(data)])/N_sim))
}

#' Estimate smoothed incidence functions and inspect deviations in the registry data.
#'
#' The first plot shows the cumulative number of diagnoses in the registry by day in black.
#' This is compared to the red line, which indicates what cumulative diagnosis would look
#' like if the rate of diagnosis was constant. The green line is a smooth fitted to the
#' actual cumulative diagnosis data. The second plot shows the deviation of the raw data from
#' the fitted smooth by day.
#'
#' The third plot shows the incidence rate per year of the registry, plotted in red. Mean
#' incidence rate is shown as a dashed black line, with the 95\% confidence interval shown with
#' dashed blue lines. A smooth fitted to the incidence data is shown in green. In the final
#' plot, the red line is the deviation of cumulative diagnoses from the fitted smooth (the same
#' as the black line in plot 2). In order to ascertain if the deviation is within reasonable
#' bounds, simulation is used. The grey lines represent N draws from a uniform distribution
#' between 0 and the last day of diagnosis, where N is the number of incident cases, plotted
#' as deviations from the smooth of cumulative diagnoses. The blue lines are 95\% confidence
#' intervals drawn from the simulated data.
#'
#' @param data A registry dataset of patient cases generated using load_data().
#' @param registry_years A vector of dates delineating years of the registry.
#' @param registry_start_year Ordinal defining the first year of the registry data to be used.
#' @param registry_end_year Ordinal defining the last year of the registry data to be used.
#' @param N Number of draws from a homogeneous Poisson process.
#' @param df Degrees of freedom for the smoothening function.
#' @return Plots of the smoothed incidence function and corresponding deviations.
#' @examples
#' smoothed_incidence(load_data(registry_data), registry_years = registry_years,
#'          registry_start_year = registry_start_year, registry_end_year = registry_end_year)
smoothed_incidence <- function(entry_date, start_date, num_years, N=1000, df=6){

  raw_incidence <- incidence(entry_date, start=start_date, num_years=num_years)

  dg <- as.numeric(difftime(entry_date, min(entry_date), units='days'))

  dfr_diags <- sort(dg)
  cum_inc <- 1:length(dfr_diags)
  smo <- smooth.spline(dfr_diags, cum_inc, df=df)

  plt1 <- plot(dfr_diags, cum_inc, pch=20, cex=0.7, xlab="days", ylab="cumulative diagnoses")
  abline(a=0, b=length(dfr_diags)/dfr_diags[length(dfr_diags)], col="red", lwd=2)
  lines(smo, col="green", lwd=2)

  plt2 <- plot(dfr_diags, cum_inc - predict(smo, dfr_diags)$y, type="l", xlab="days", ylab="deviation from smooth")
  mean_rate <- mean(raw_incidence)
  day_mean_rate <- mean_rate/365
  CI_lim <- 1.96 * sqrt(mean_rate)/365
  pl_lim <- CI_lim * 2.0
  pre_smo <- predict(smo, 1:(365*num_years), deriv=1)

  plt3 <- plot(365*(1:num_years) - 182.5, raw_incidence/365, pch=20, col="red",
       xlab="days", ylab="incidence rate",
       ylim=c(day_mean_rate-pl_lim, day_mean_rate+pl_lim ))
  lines(365*(1:num_years) - 182.5, raw_incidence/365, col="red",lwd=2)
  lines(pre_smo, type="l", lwd=2, col="green")

  abline(h = day_mean_rate, lty=2)
  abline(h = day_mean_rate - CI_lim, lty=3, col="blue")
  abline(h = day_mean_rate + CI_lim, lty=3, col="blue")
  abline(v=(1:num_years)*365, col="pink", lty=2)

  N <- length(dfr_diags)
  M <- 1000
  boot_out <- matrix(NA, nrow = M, ncol = N)

  for (i in 1:M){
    x <- sort(runif(N, 0, max(dfr_diags)))
    the_smo <- smooth.spline(x, 1:N, df=4)
    boot_out[i, ] <- (1:N) - predict(the_smo, x)$y
  }

  plt4 <- plot(NA, xlim=c(0,max(dfr_diags)), ylim=c(-20,20), xlab="days", ylab="deviation from smooth")
  for (i in 1:M){
    lines(x, boot_out[i,], col="grey")
  }

  lines(dfr_diags, cum_inc - predict(smo, dfr_diags)$y, col="red")

  upper_lim <- apply(boot_out, 2, quantile, probs=0.975)
  lower_lim <- apply(boot_out, 2, quantile, probs=0.025)

  lines(x, upper_lim, col="blue")
  lines(x, lower_lim, col="blue")

  return(list(plt1, plt2, plt3, plt4))

}


smoothed_incidence_current <- function(entry_date, start_date, num_years, N=1000, df=6){

  raw_incidence <- incidence(entry_date, start=start_date, num_years=num_years)

  dg <- as.numeric(difftime(entry_date, min(entry_date), units='days'))

  dfr_diags <- sort(dg)
  cum_inc <- 1:length(dfr_diags)
  smo <- smooth.spline(dfr_diags, cum_inc, df=df)

  plt1 <- plot(dfr_diags, cum_inc, pch=20, cex=0.7, xlab="days", ylab="cumulative diagnoses")
  abline(a=0, b=length(dfr_diags)/dfr_diags[length(dfr_diags)], col="red", lwd=2)
  lines(smo, col="green", lwd=2)

  plt2 <- plot(dfr_diags, cum_inc - predict(smo, dfr_diags)$y, type="l", xlab="days", ylab="deviation from smooth")
  mean_rate <- mean(raw_incidence)
  day_mean_rate <- mean_rate/365
  CI_lim <- 1.96 * sqrt(mean_rate)/365
  pl_lim <- CI_lim * 2.0
  pre_smo <- predict(smo, 1:(365*num_years), deriv=1)

  plt3 <- plot(365*(1:num_years) - 182.5, raw_incidence/365, pch=20, col="red",
       xlab="days", ylab="incidence rate",
       ylim=c(day_mean_rate-pl_lim, day_mean_rate+pl_lim ))
  lines(365*(1:num_years) - 182.5, raw_incidence/365, col="red",lwd=2)
  lines(pre_smo, type="l", lwd=2, col="green")

  abline(h = day_mean_rate, lty=2)
  abline(h = day_mean_rate - CI_lim, lty=3, col="blue")
  abline(h = day_mean_rate + CI_lim, lty=3, col="blue")
  abline(v=(1:num_years)*365, col="pink", lty=2)

  N <- length(dfr_diags)
  M <- 1000
  boot_out <- matrix(NA, nrow = M, ncol = N)

  for (i in 1:M){
    x <- sort(runif(N, 0, max(dfr_diags)))
    the_smo <- smooth.spline(x, 1:N, df=4)
    boot_out[i, ] <- (1:N) - predict(the_smo, x)$y
  }

  plt4 <- plot(NA, xlim=c(0,max(dfr_diags)), ylim=c(-20,20), xlab="days", ylab="deviation from smooth")
  for (i in 1:M){
    lines(x, boot_out[i,], col="grey")
  }

  lines(dfr_diags, cum_inc - predict(smo, dfr_diags)$y, col="red")

  upper_lim <- apply(boot_out, 2, quantile, probs=0.975)
  lower_lim <- apply(boot_out, 2, quantile, probs=0.025)

  lines(x, upper_lim, col="blue")
  lines(x, lower_lim, col="blue")

  return(list(plt1, plt2, plt3, plt4))

}

smoothed_incidence_gg <- function(entry_date, start_date, num_years, N=1000, df=6){

  raw_incidence <- incidence(entry_date, start=start_date, num_years=num_years)

  dg <- as.numeric(difftime(entry_date, min(entry_date), units='days'))

  dfr_diags <- sort(dg)
  cum_inc <- 1:length(dfr_diags)
  smo <- smooth.spline(dfr_diags, cum_inc, df=df)

  plt1 <- plot(dfr_diags, cum_inc, pch=20, cex=0.7, xlab="days", ylab="cumulative diagnoses")
  abline(a=0, b=length(dfr_diags)/dfr_diags[length(dfr_diags)], col="red", lwd=2)
  lines(smo, col="green", lwd=2)

  plt2 <- plot(dfr_diags, cum_inc - predict(smo, dfr_diags)$y, type="l", xlab="days", ylab="deviation from smooth")
  mean_rate <- mean(raw_incidence)
  day_mean_rate <- mean_rate/365
  CI_lim <- 1.96 * sqrt(mean_rate)/365
  pl_lim <- CI_lim * 2.0
  pre_smo <- predict(smo, 1:(365*num_years), deriv=1)

  plt3 <- plot(365*(1:num_years) - 182.5, raw_incidence/365, pch=20, col="red",
       xlab="days", ylab="incidence rate",
       ylim=c(day_mean_rate-pl_lim, day_mean_rate+pl_lim ))
  lines(365*(1:num_years) - 182.5, raw_incidence/365, col="red",lwd=2)
  lines(pre_smo, type="l", lwd=2, col="green")

  abline(h = day_mean_rate, lty=2)
  abline(h = day_mean_rate - CI_lim, lty=3, col="blue")
  abline(h = day_mean_rate + CI_lim, lty=3, col="blue")
  abline(v=(1:num_years)*365, col="pink", lty=2)

  N <- length(dfr_diags)
  M <- 1000
  boot_out <- matrix(NA, nrow = M, ncol = N)

  for (i in 1:M){
    x <- sort(runif(N, 0, max(dfr_diags)))
    the_smo <- smooth.spline(x, 1:N, df=4)
    boot_out[i, ] <- (1:N) - predict(the_smo, x)$y
  }

  plt4 <- plot(NA, xlim=c(0,max(dfr_diags)), ylim=c(-20,20), xlab="days", ylab="deviation from smooth")
  for (i in 1:M){
    lines(x, boot_out[i,], col="grey")
  }

  lines(dfr_diags, cum_inc - predict(smo, dfr_diags)$y, col="red")

  upper_lim <- apply(boot_out, 2, quantile, probs=0.975)
  lower_lim <- apply(boot_out, 2, quantile, probs=0.025)

  lines(x, upper_lim, col="blue")
  lines(x, lower_lim, col="blue")

  return(list(plt1, plt2, plt3, plt4))

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

#' Plot the age distribution of incident cases in the registry data.
#'
#' @param data A registry dataset of patient cases generated using load_data().
#' @param registry_years A vector of dates delineating years of the registry.
#' @param registry_start_year Ordinal defining the first year of the registry data to be used.
#' @param registry_end_year Ordinal defining the last year of the registry data to be used.
#' @param df Degrees of freedom for the smooth.
#' @return Plot of the raw data and smoothed age distribution function.
#' @examples
#' incidence_age_distribution_current(load_data(registry_data), registry_years, registry_start_year = registry_start_year,
#' registry_end_year = registry_end_year)
incidence_age_distribution_current <- function(data, registry_years, registry_start_year, registry_end_year, df=10){
    
    ages <- rep(0, 100)
    for (i in 1:dim(data)[1]){
        ages[1 + floor(data$age_initial[i])] <- ages[1 + floor(data$age_initial[i])] + 1
    }
    plot(0:99, ages[1:100], pch=20, xlab="age (years)", ylab="incident cases per year")
    smage <- smooth.spline(0:99, ages[1:100], df=df)
    lines(smage, col="blue", lwd=2)
    
}

#' Inspect consistency of survival data between years of the registry and with a Cox Proportional Hazards model.
#'
#' The first plot is of the Kaplan-Meier survival curve on total cases in the registry. The second plot is the
#' Kaplan-Meier survival curve for each age group, as delineated by the user using the "ages" argument. The third
#' plot is of the residuals between the raw data and the fitted Cox Proportional Hazards model. If the model is
#' a good representation of the data the line should be horizontal. The last plot is of Kaplain Meier survival
#' curves fitted to cases subdivided by year of diagnosis within the registry (black lines), compared to total
#' cases shown in blue.
#'
#' @param data A registry dataset of patient cases generated using load_data().
#' @param ages A vector of ages at which to break the dataset for Kaplan Meier plotting.
#' @param registry_years A vector of dates delineating years of the registry.
#' @param registry_start_year Ordinal defining the first year of the registry data to be used.
#' @param registry_end_year Ordinal defining the last year of the registry data to be used.
#' @return A sequence of plots indicated the consistency of survival data between years of the registry and with a Cox Proportional Hazards model.
#' @examples
#' survival_modelling_diagnostics(load_data(registry_data), registry_years, registry_start_year = registry_start_year, registry_end_year = registry_end_year, ages = c(55, 65, 75, 85, 100))
survival_modelling_diagnostics <- function(data, ages, registry_years, registry_start_year,
                      registry_end_year){

  if (is.numeric(ages) != TRUE) stop("error: ages is not numeric.")
  if (is.vector(ages) != TRUE) stop("error: ages is not a vector.")

  years_estimated <- registry_end_year - registry_start_year + 1

  dfr_r <- data[data$date_initial >= registry_years[registry_start_year], ]

  km <- survfit(Surv(survival_time, indicator) ~ 1, data=dfr_r)
  plt1 <- plot(km, lwd=2, col="blue", xlab="survival (days)", ylab="probability")

  km2 <- survfit(Surv(survival_time, indicator) ~ cut(age_initial, breaks=ages), data=dfr_r)
  plt2 <- plot(km2, lwd=2, col=1:length(ages), xlab="survival (days)", ylab="probability")

  cx <- coxph(Surv(survival_time, indicator) ~ age_initial, data=dfr_r)

  halves <- rep(0, length(ages) - 1)
  for(i in 1:length(halves)){
    halves[i] <- mean(c(ages[i], ages[i + 1]))
  }

  cxp <- survfit(cx, newdata=data.frame(age_initial=halves))
  lines(cxp, lwd=2, col=1:6, lty=2, mark.time=F)

  output <- cox.zph(cx)
  plt3 <- plot(cox.zph(cx))

  plt4 <- plot(km, lwd=2, col="blue", mark.time=F, conf.int=T, xlab="survival (days)", ylab="probability")

  for (i in registry_start_year:registry_end_year){
    dfr_L <- dfr_r[dfr_r$date_initial >= registry_years[i] & dfr_r$date_initial < registry_years[i + 1], ]
    kmlines <- survfit(Surv(survival_time, indicator) ~ 1, data=dfr_L)
    lines(kmlines, mark.time=F, conf.int=F)
  }

  return(list(plt1, plt2, plt3, plt4, output))

}

#' Inspect functional form of age.
#'
#' @param data A registry dataset of patient cases.
#' @param df Degrees of freedom for the smooth.
#' @return Plots of the functional form of age.
#' @examples
#' functional_form_age(registry_data_r)
functional_form_age <- function(data, df=4){
    
  ### TO DO
  # Parse formula for survival object and tweak all the variable names
  # Figure out how to do a good unit test with the previous version - setting seed before calling coxph didn't work
  ###
    
  set.seed(17)
  cxnl <- coxph(Surv(survival_time, indicator) ~ pspline(age_initial, df=df), data)
  output1 <- summary(cxnl)

  plt1 <- termplot(cxnl)

  f <<- datadist(data)
  options(datadist="f")

  mod_rms <- cph(Surv(survival_time, indicator) ~ rcs(age_initial, df), data, x=TRUE, y=TRUE, surv=T, time.inc=1)
  output2 <- anova(mod_rms)
  output3 <- summary(mod_rms)

  plt2 <- plot(Predict(mod_rms, age_initial), lwd=3, adj.subtitle=T)

  return(list(plt1, plt2, output1, output2, output3))

}

#' Inspect functional form of age.
#'
#' @param data A registry dataset of patient cases generated using load_data().
#' @param registry_years A vector of dates delineating years of the registry.
#' @param registry_start_year Ordinal defining the first year of the registry data to be used.
#' @param registry_end_year Ordinal defining the last year of the registry data to be used.
#' @param df Degrees of freedom for the smooth.
#' @return Plots of the functional form of age.
#' @examples
#' functional_form_age_current(load_data(registry_data), registry_years, registry_start_year = registry_start_year, registry_end_year = registry_end_year)
functional_form_age_current <- function(data, registry_years, registry_start_year, registry_end_year, df=4){
    
    data_r <- data[data$date_initial >= registry_years[registry_start_year] & data$date_initial < registry_years[registry_end_year], ]
    set.seed(17)    
    cxnl <- coxph(Surv(survival_time, indicator) ~ pspline(age_initial, df=df), data=data_r)
    output1 <- summary(cxnl)
    
    plt1 <- termplot(cxnl)
    
    f <<- datadist(data_r)
    options(datadist="f")
    
    mod_rms <- cph(Surv(survival_time, indicator) ~ rcs(age_initial, df), data=data_r, x=TRUE, y=TRUE, surv=T, time.inc=1)
    output2 <- anova(mod_rms)
    output3 <- summary(mod_rms)
    
    plt2 <- plot(Predict(mod_rms, age_initial), lwd=3, adj.subtitle=T)
    
    return(list(plt1, plt2, output1, output2, output3))
    
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
#' @param num_years Integer representing the number of complete years of the registry for which incidence is to be calculated.
#' @return Chi-squared test of difference between prevalence prediction and counted prevalence at the index date.
#' @examples
#' prev_chisq(entry = registry_data$entrydate,
#'            events = registry_data$eventdate, 
#'            status = registry_data$status, 
#'            start="2004-01-30", num_years = 9,
#'            by_year = by_year_total)
prev_chisq <- function(entry, events, status, by_year, start=NULL, num_years=NULL){
  observed <- counted_prevalence(entry, events, status, start, num_years)
  predicted <- rev(by_year[1:num_years])
  chi <- sum(((observed - predicted)^2)/predicted)
  1 - pchisq(chi, num_years - 1)

}

#' Chi squared test between prevalence prediction and observed values in the registry.
#'
#' @param data A registry dataset of patient cases generated using load_data().
#' @param registry_years A vector of dates delineating years of the registry.
#' @param registry_start_year Ordinal defining the first year of the registry data to be used.
#' @param registry_end_year Ordinal defining the last year of the registry data to be used.
#' @return Plots of the smoothed incidence function and corresponding deviations.
#' @examples
#' prev_chisq_current(load_data(registry_data), registry_years = registry_years, registry_start_year = registry_start_year, registry_end_year = registry_end_year, by_year = by_year_total)
prev_chisq_current <- function(data, registry_years, registry_start_year, registry_end_year, by_year){
    
    observed <- counted_prevalence_current(data, registry_years, registry_start_year, registry_end_year)
    predicted <- rev(by_year[1:(registry_end_year - registry_start_year + 1)])
    chi <- sum(((observed - predicted)^2)/predicted)
    result <- 1 - pchisq(chi, (registry_end_year - registry_start_year))
    return(result)
    
}