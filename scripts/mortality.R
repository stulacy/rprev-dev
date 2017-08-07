# Interpolates daily mortality rates for every age and sex combination from the
# yearly values in UKmortality

library(dplyr)
library(rprev)
data(UKmortality)

interpolate_daily_mortality <- function(data) {
    DAYS_IN_YEAR <- 365.25

    # Could probably improve on this extraction of response and variable
    rate <- vapply(seq(0, max_age-1),
                   function(x) mean(data[floor(data$age) == x, 'rate']),
                   numeric(1))
    num_ages <- length(rate)

    a_rate <- c(2 * rate[1] - rate[2], rate, 2 * rate[num_ages] - rate[num_ages - 1])
    base <- 183 + DAYS_IN_YEAR * (0:num_ages) # Where does 183/182 come from?
    base <- c(-182, base)
    daily_rate <- approx(base, a_rate, 1:(num_ages * DAYS_IN_YEAR))
    daily_rate$y / DAYS_IN_YEAR
}

vals <- data.frame(expand.grid(age=0:100, sex=c(0,1)))
daily_rate <- data.frame(setNames(lapply(c(0, 1), function(s) {
    df_sub <- filter(UKmortality, sex==s)
    interpolate_daily_mortality(df_sub)
}), c('M', 'F'))) %>%
    mutate(overall=(M+F)/2,
           agedays = 1:nrow(.))

save(daily_rate, file="data/UKmortalitydays.rda")
