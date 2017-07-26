library(flexsurvcure)
library(rprev)
library(ggplot2)
library(tidyr)
data(prevsim)

cure_model <- flexsurvcure(Surv(time, status) ~ age + sex, data=prevsim, link="logistic",
                           dist='weibull', mixture=T)

# I assume theta is the baseline probability of being susceptible?
cure_model

# Does this mean that covariates aren't being placed on the distribution location or scale parameters then?
# The docs say can use anc to specify effects on ancilliary covariates but does this just mean the scale parameters rather than the location?
# Ideally want to know exactly what is being modelled here


### CALCULATING CDF

# See my paper notes for efforts to determine inverse of CDF for drawing event times
# Let's plot the S(t) and thereby CDF for an individual with known attributes.
# Since the flexsurvcure package produces relative survival estimates we'll simulate population mortality

# Simulate 11 years of population mortality
pop_mort_raw <- cumsum(runif(11))
pop_mort_norm <- 1 - (pop_mort_raw / max(pop_mort_raw)) * 0.3
pop_mort_norm

# Expand into daily values
daily_pop_mort <- rep(pop_mort_norm, each=365)
length(daily_pop_mort)

# Now predict relative survival for 4016 days for a 50 year old man
pred <- summary(cure_model, newdata=data.frame(age=50, sex=factor(1)), t=1:4015)
pred_surv <- pred$`age=50, sex=1`$est
pred_surv_overall <- pred_surv * daily_pop_mort

data.frame(time=1:length(pred_surv),
           rel = pred_surv,
           overall = pred_surv_overall) %>%
    gather(method, prob, -time) %>%
    ggplot(aes(x=time, y=prob, colour=method)) +
        geom_line() +
        ylim(0, 1) +
        theme_bw()

# Yep that is 100% a discrete survival curve...

# So how do we draw from it?!


