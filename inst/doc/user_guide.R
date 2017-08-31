## ----setup, message = FALSE, warning = FALSE-----------------------------
library(rprev)
library(survival)
data(prevsim)

## ------------------------------------------------------------------------
summary(prevsim)

## ---- echo=F-------------------------------------------------------------
#As our registry is short relative to the duration of survival from this disease, we do not have much information about long-term survival after this point in our patient population in order to base a model. In this scenario, instead of extrapolating survival (and thus risking the problem of "immortal" subjects), we recommend using general population survival data to model long-term survivors, and we refer to this as a *cure model* [@crouch2014determining]. In this model, after a specified "cure time", surviving subjects are considered to revert to the survival characteristics for someone of their age and sex in the general population. In our example, a patient is considered *cured* after five years survival with the disease.

#Until the point in time that the cure model takes effect, survival is modelled on disease registry data; after this time, survival probability for the surviving cases is modelled using population mortality rates. `rprev` includes data from the Office of National Statistics (ONS) [see @ons, data licenced under the Open Government Licence v3.0] to describe UK mortality rates (in data set `UKmortality`), however, the user may supply their own population mortality data. 

## ----basicsurvival, fig.height=4, fig.width=7, echo=F--------------------
survf <- survfit(Surv(time, status) ~ sex, data=prevsim)
survf_df <- data.frame(t=survf$time, s=survf$surv, sex=rep(c('M', 'F'), survf$strata))
ggplot2::ggplot(survf_df, aes(x=t, y=s, colour=sex)) +
    ggplot2::geom_line() +
    theme_bw() +
    labs(x='Time (days)', y='Survival probability') +
    ylim(0, 1)

## ------------------------------------------------------------------------
inc <- incidence(prevsim$entrydate, population_size=3.2e6, start='2005-09-01', num_reg_years=5)

## ------------------------------------------------------------------------
inc
names(inc)

## ------------------------------------------------------------------------
summary(inc)

## ----incidencerate, error = TRUE-----------------------------------------
inc$mean

## ---- fig.width = 7, fig.height = 4--------------------------------------
plot(inc)

## ----incidenceage, fig.width = 7, fig.height = 4, error = TRUE-----------
prevsim_r <- prevsim[prevsim$entrydate >= "2004-01-30", ]

incidence_age_distribution(prevsim_r$age)

## ----survivaldiag, fig.width = 7, fig.height = 4-------------------------
km <- survfit(Surv(time, status) ~ 1, data=prevsim_r)
plot(km, lwd=2, col="blue", main="Overall Survival", xlab="Days", 
     ylab="Survival probability")

## ----survivaldiag2, fig.width = 7, fig.height = 4------------------------
ages = c(55, 65, 75, 85, 100)
km2 <- survfit(Surv(time, status) ~ cut(age, breaks=ages), data=prevsim_r)
plot(km2, lwd=2, col=1:length(ages), main="Survival stratified by age", xlab="Days", 
     ylab="Survival probability")
legend("topright", legend=substring(names(km2$strata), 25, 32), lty = 1, 
       col=1:length(ages))

## ----survivaldiag4, fig.width = 7, fig.height = 4, results='hide'--------
plot(km, lwd=2, col="blue", mark.time=F, conf.int=T, xlab="Days", 
     ylab="Survival probability")
num_reg_years <- 9
registry_years <- determine_yearly_endpoints(date='2004-01-30', 
                                             num_years=num_reg_years)
sapply(seq(num_reg_years),
       function(i) lines(survfit(Surv(time, status) ~ 1, 
                                 data=prevsim_r[prevsim_r$entrydate >= 
                                                          registry_years[i] & 
                                                          prevsim_r$entrydate < 
                                                          registry_years[i + 1], ]), 
                         mark.time = F, conf.int = F))

## ----survivaldiag3, fig.width = 7, fig.height = 4------------------------
cx <- coxph(Surv(time, status) ~ age, data=prevsim_r)
cxp <- survfit(cx, 
               newdata=data.frame(age=sapply(seq(length(ages) - 1), 
                                             function(i) mean(c(ages[i], ages[i + 1]))))) 
plot(cox.zph(cx))
lines(cxp, lwd=2, col=1:length(ages), lty=2, mark.time=F)


## ------------------------------------------------------------------------
cox.zph(cx)

## ----ageform2, fig.width = 7, fig.height = 4, error = TRUE---------------
functional_form_age(Surv(time, status) ~ age, prevsim_r, df=4, plot_fit=T)

## ----weibull, error = TRUE-----------------------------------------------
wb <- survreg(Surv(time, status) ~ age + sex, data=prevsim_r)
summary(wb)

## ----popsurv, fig.width = 7, fig.height = 4, error = TRUE----------------
data(UKmortality)

daily_survival_males <- population_survival_rate(rate ~ age, subset(UKmortality, 
                                                                    sex == 'M'))
daily_survival_females <- population_survival_rate(rate ~ age, subset(UKmortality, 
                                                                      sex == 'F'))

plot(daily_survival_males, type="l", col="blue", xlab="days", ylab="survival")
lines(daily_survival_females, col="red")
legend("topright", legend = c("Males", "Females"),
               bty = "n", lty = 1, col = c(4,2))

## ----survival, error = TRUE----------------------------------------------
prevalence_counted_legacy(prevsim$entrydate, 
                          prevsim$eventdate, 
                          prevsim$status, 
                          index_date="2013-01-30",
                          num_reg_years=9) 

## ----prevalencetotal, error = TRUE---------------------------------------
prevalence_total <- prevalence_legacy(Surv(time, status) ~ sex(sex) + age(age) + 
                                   entry(entrydate) + event(eventdate),
                               prevsim, num_years_to_estimate=c(3, 5, 10), 
                               index_date='2013-01-30',
                               num_reg_years=9, cure=5) 

## ------------------------------------------------------------------------
prevalence_total

## ------------------------------------------------------------------------
summary(prevalence_total)

## ------------------------------------------------------------------------
prevalence_total$estimates

## ------------------------------------------------------------------------
prevsurv <- survfit(prevalence_total, newdata=list(age=60, sex=0))
prevsurv

## ------------------------------------------------------------------------
summary(prevsurv, years=c(1, 3, 5, 10))

## ---- fig.width=7, fig.height=4------------------------------------------
plot(prevsurv, pct_show=0.90)

## ------------------------------------------------------------------------
prevalence_total$pval

## ----test, error = TRUE--------------------------------------------------
test_prevalence_fit(prevalence_total)

## ---- fig.height = 4, fig.width = 7--------------------------------------
hist(prevsim$age[prevsim$entrydate >= min(registry_years) & 
                     prevsim$entrydate < max(registry_years)], 
     col=rgb(1,0,0, alpha=0.5), xlim=c(0,100), ylim=c(0,0.045), 
     main = "", xlab = "age", prob = TRUE)
hist(unlist(prevalence_total$simulated$posterior_age),
     col=rgb(0,1,0, alpha=0.5), prob = TRUE, add=T)
legend("topleft", legend = c("Incident", "Prevalent"), bty = "n", lty = 1, 
       col = c(rgb(1,0,0, alpha=0.5),
               rgb(0,1,0, alpha=0.5)))

