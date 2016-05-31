#' ...
#'
#' @param age Age of example patient.
#' @param sex Sex of example patient ("Male" or "Female").
#' @param wb Weibull model fitted to the data.
#' @return ...
.dfr_sc <- function(age, sex, wb) {
    if (sex == "Male") {
        sexn <- 0
    } else {
        sexn <- 1
    }
    return(exp(wb$coef[1] + age*wb$coef[2] + sexn*wb$coef[3]))
}

#' ...
#'
#' @param t ...
#' @param shape ...
#' @param scale ...
#' @return ...
.wb_Su <- function(t, shape, scale) {
    return(exp(-(t/scale)^shape))
}

#' ...
#'
#' @param t ...
#' @param shape ...
#' @param scale ...
#' @return ...
.wb_haz <- function(t, shape, scale) {
    return((shape/scale)*(t/scale)^(shape-1))
}

#' ...
#'
#' @param age Age of example patient.
#' @param sex Sex of example patient ("Male" or "Female").
#' @param the_length Number of days to predict.
#' @param wb Weibull model fitted to the data.
#' @return ...
.plot_haz <- function(age, sex, the_length=3000, wb) {
    s_time <- seq(0, the_length, by=1)
    haz <- wb_haz(s_time, dfr_sh, dfr_sc(age, sex, wb))
    plot(s_time, haz, type="l", lwd=2, col="red", ylim=c(0,0.002),
         main=paste("age = ", age, ", sex =", sex), xlab="survival (days)", ylab="probability")
}

#' ...
#'
#' @param data Population survival dataset.
#' @param age Age of example patient.
#' @param sex Sex of example patient ("Male" or "Female").
#' @param the_length Number of days to predict.
#' @return ...
.plot_population_hazard <- function(data, age, sex, the_length=3000) {
    if(sex == "Male") sex = "Males"
    if(sex == "Female") sex = "Females"
    daily_survival_rate <- function(data, sex)
        lines(0:(the_length-1), daily_survival_rate[(365*age):(365*age + the_length - 1)],
              col="blue", lwd=2)
}

#' ...
#'
#' @param age Age of example patient.
#' @param sex Sex of example patient ("Male" or "Female").
#' @param the_length Number of days to predict.
#' @param wb Weibull model fitted to the data.
#' @return ...
.plot_Su <- function(age, sex, the_length=3000, wb) {
    s_time <- seq(0, the_length, by=1)
    Su <- wb_Su(s_time, dfr_sh, dfr_sc(age, sex, wb))
    plot(s_time, Su, type="l", lwd=2, col="red", ylim=c(0,1),
         main=paste("age = ", age, ", sex =", sex), xlab="survival (days)", ylab="probability")
}

#' ...
#'
#' @param data Population survival dataset.
#' @param age Age of example patient.
#' @param sex Sex of example patient ("Male" or "Female").
#' @param the_length Number of days to predict.
#' @return ...
.plot_pop_Su <- function(data, age, sex, the_length=3000) {
    if(sex == "Male") sex = "Males"
    if(sex == "Female") sex = "Females"
    daily_survival_rate <- function(data, sex)
        lines(1:the_length, cumprod(1 - daily_survival_rate[(age*365):(age*365 + the_length - 1)]),
              col="cyan", lwd=2)
}

#' ...
#'
#' @param data A registry dataset of patient cases generated using load_data().
#' @param registry_years A vector of dates delineating years of the registry.
#' @param registry_start_year Ordinal defining the first year of the registry data to be used.
#' @param age Age of example patient.
#' @param sex Sex of example patient ("Male" or "Female").
#' @param limits ...
#' @param colour ...
#' @return ...
.plot_km <- function(data, registry_years, registry_start_year, age,
                    sex, limits, colour="green") {
    if (sex == "Male") {
        sexn <- 0
    } else {
        sexn <- 1
    }
    dfr_r <- data[data$date_initial >= registry_years[registry_start_year],]
    kma <- survival::survfit(Surv(survival_time, indicator) ~ 1,
                             data=dfr_r[dfr_r$sex == sexn & dfr_r$age_initial > age-limits &
                             dfr_r$age_initial < age + limits, ])
    lines(kma, lwd=1, col=colour, conf.int=T)
}