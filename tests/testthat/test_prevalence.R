library(prevR)
library(survival)
library(abind)
context('Prevalence')
data(prevsim)

test_that("prevalence returns same values as before", {
    set.seed(3)
    i <- 1
    expect_ref <- function(data, num_years_to_estimate, start, years, cure, boot) {
        fn <- paste('cache/prevalence/prevalence_', i, '.rds', sep='')
        expect_equal_to_reference(prevalence(Surv(time, status) ~ sex(sex) + age(age) + entry(entrydate),
                                             data, num_years_to_estimate,
                                             start=start, num_reg_years=years,
                                             cure=cure, N_boot=boot), file=fn)
        i <<- i + 1
    }
    expect_ref(prevsim, 10, '2004-01-01', 9, cure=5, boot=10)
    expect_ref(prevsim, 10, '2004-01-01', 9, cure=5, boot=20)
    expect_ref(prevsim, 10, '2005-05-23', 4, cure=5, boot=10)
    expect_ref(prevsim, 10, '2004-01-01', 9, cure=3, boot=10)
    expect_ref(prevsim, 4, '2004-01-01', 3, cure=5, boot=10)
    expect_ref(prevsim, 4, '2004-01-01', 3, cure=1, boot=10)
})

test_that("prevalence with a thousand bootstraps returns same values as before", {
    skip_on_cran()
    skip("too slow")
    set.seed(3)
    i <- 1
    expect_ref <- function(data, num_years_to_estimate, start, years, cure, boot) {
        fn <- paste('cache/prevalence/prevalence_thousand_', i, '.rds', sep='')
        expect_equal_to_reference(prevalence(Surv(time, status) ~ sex(sex) + age(age) + entry(entrydate),
                                             data, num_years_to_estimate,
                                             start=start, num_reg_years=years,
                                             cure=cure, N_boot=boot), file=fn)
        i <<- i + 1
    }
    expect_ref(prevsim, 10, '2004-01-01', 9, cure=5, boot=1000)
    expect_ref(prevsim, 10, '2006-09-03', 5, cure=5, boot=1000)
})

# TODO error testing, that if num_years_to_estimate is longer than a set amount it throws an error

test_that("Error is raised when passing a population data frame not setup correctly", {
    expect_poperror <- function(popdata, msg) {
        expect_error(prevalence(Surv(time, status) ~ sex(sex) + age(age) + entry(entrydate),
                                             prevsim, num_years_to_estimate=10, population_data=popdata,
                                             start='2004-01-01', num_reg_years=9,
                                cure=5, N_boot=10), msg)
    }
    missing_cols <- "Error: the supplied population data frame must contain columns 'rate', 'age', 'sex'."
    missing_levels <- "Error: the same levels must be present in both the population and registry data. '0' and '1' by default where male is 0."

    expect_poperror(data.frame(age=rnorm(10, 60, 10), gender=rbinom(10, 1, 0.5)), missing_cols)  # Lacks rate and sex is named as gender
    expect_poperror(data.frame(age=rnorm(10, 60, 10), gender=rbinom(10, 1, 0.5), rate=runif(10)), missing_cols)  # sex is named as gender
    expect_poperror(data.frame(gender=rbinom(10, 1, 0.5), rate=runif(10)), missing_cols)  # Lacking age and sex is misnamed
    expect_poperror(data.frame(sex=rbinom(10, 1, 0.5), rate=runif(10)), missing_cols)  # Lacking age
    expect_poperror(data.frame(sex=rbinom(10, 1, 0.5)), missing_cols)  # Lacking age and rate
    expect_poperror(data.frame(age=rnorm(10, 60, 10), sex=sample(c('M', 'F'), 10, replace=T), rate=runif(10)), missing_levels)  # sex has different levels to that in prevsim
})

test_that("Formula for prevalence must have age, sex, and entry functions.", {
    expect_formerror <- function(form, msg) {
        expect_error(prevalence(form,
                                prevsim, num_years_to_estimate=10,
                                start='2004-01-01', num_reg_years=9,
                                cure=5, N_boot=10), msg)
    }
    missing_funcs <- "Error: provide function terms for age, sex, and entry date."
    expect_formerror(Surv(time, status) ~ age + sex + entry, missing_funcs)
    expect_formerror(Surv(time, status) ~ age(age) + sex + entry, missing_funcs)
    expect_formerror(Surv(time, status) ~ age + sex(sex) + entr, missing_funcs)
    expect_formerror(Surv(time, status) ~ age + sex + entry(entry), missing_funcs)
    expect_formerror(Surv(time, status) ~ age(age) + sex + entry(entry), missing_funcs)

    # All the functions are present but the columns aren't! i.e. entry is entrydate in prevsim dataset
    missing_cols <- "undefined columns selected"
    expect_formerror(Surv(time, status) ~ age(age) + sex(sex) + entry(entry), missing_cols)
    expect_formerror(Surv(time, status) ~ age(age) + sex(gender) + entry(entrydate), missing_cols)

})

test_that("Error is raised when levels for sex aren't the same in registry and population data", {
    expect_sexerror <- function(regdata, popdata, msg) {
        expect_error(prevalence(Surv(time, status) ~ sex(sex) + age(age) + entry(entrydate),
                                             regdata, num_years_to_estimate=10, population_data=popdata,
                                             start='2004-01-01', num_reg_years=9,
                                cure=5, N_boot=10), msg)
    }
    missing_levels <- "Error: the same levels must be present in both the population and registry data. '0' and '1' by default where male is 0."

    regN = 1000
    popN = 5000
    agereg <- rnorm(regN, 60, 10)
    agepop <- rnorm(popN, 55, 20)
    rate <- runif(popN)
    time <- runif(regN, 0, 1000)
    status <- rbinom(regN, 1, 0.8)
    entry <- as.character(sapply(sample(seq(8), regN, replace=T), function(x) sprintf("200%d-03-17", x)))

    pop_01 <- data.frame(age=agepop, sex=rbinom(10, 1, 0.5), rate=rate)
    pop_mf <- data.frame(age=agepop, sex=sample(c('M', 'F'), 10, replace=T), rate=rate)

    reg_01 <- data.frame(time=time, status=status, age=agereg, sex=rbinom(10, 1, 0.5), entrydate=entry)
    reg_mf <- data.frame(time=time, status=status, age=agereg, sex=sample(c('M', 'F'), 10, replace=T), entrydate=entry)

    expect_sexerror(reg_01, pop_mf, missing_levels)
    expect_sexerror(reg_mf, pop_01, missing_levels)
})

test_that("Prevalence messages the user when the number of registry years are greater than the number of years asked to predict prevalence for", {
    expect_msg <- function(Nyears, start, regyears=NULL) {
        expect_message(prevalence(Surv(time, status) ~ sex(sex) + age(age) + entry(entrydate),
                                  prevsim, num_years_to_estimate=Nyears,
                                  start=start, num_reg_years=regyears,
                                  cure=5, N_boot=10))
    }
    expect_msg(Nyears=5, start='2003-01-01', regyears=10)
    expect_msg(Nyears=5, start='2003-01-01')  # Note there are 10 registry years in this dataset which are used if not provided
    expect_msg(Nyears=1, start='2011-01-01')  # However if we set registry to start at 2011 then only 2 years of registry data
})

test_that("counted_prevalence results in the same values as before", {
    set.seed(3)
    i <- 1
    expect_ref <- function(entry, event, status, start, years=NULL, index=NULL) {
        fn <- paste('cache/prevalence/countedprevalence_', i, '.rds', sep='')
        expect_equal_to_reference(counted_prevalence(entry, event, status, start=start,
                                                     num_reg_years=years), file=fn)
        i <<- i + 1
    }
    expect_ref(prevsim$entrydate, prevsim$eventdate, prevsim$status, '2004-03-05', 5)
    expect_ref(prevsim$entrydate, prevsim$eventdate, prevsim$status, '2004-03-05')
})

test_that("n_year_estimates returns the same values as before", {
    set.seed(3)
    i <- 1
    expect_ref <- function(data, num_years_to_estimate, start, years, cure, boot, nyearest) {
        fn <- paste('cache/prevalence/nyearestimates_', i, '.rds', sep='')
        obj <- prevalence(Surv(time, status) ~ sex(sex) + age(age) + entry(entrydate),
                          data, num_years_to_estimate,
                          start=start, num_reg_years=years,
                          cure, N_boot=boot)
        expect_equal_to_reference(n_year_estimates(obj, num_years_to_estimate=nyearest,
                                                   population_size = 35e5), file=fn)
        i <<- i + 1
    }

    expect_ref(prevsim, 10, '2004-01-01', 9, cure=5, boot=10, nyearest=10)
    expect_ref(prevsim, 10, '2005-04-21', 8, cure=3, boot=10, nyearest=3)
    expect_ref(prevsim, 5, '2006-04-21', 3, cure=3, boot=10, nyearest=3)
})

test_that("n_year_estimates correctly throws an error when asked to estimate confidence intervals for more years than initially estimated.", {
    msg = "Error: can't calculate prevalence for more years than present in the prevalence object."
    expect_err <- function(data, num_years_to_estimate, start, years, cure, boot, nyearest) {
        obj <- prevalence(Surv(time, status) ~ sex(sex) + age(age) + entry(entrydate),
                          data, num_years_to_estimate,
                          start=start, num_reg_years=years,
                          cure, N_boot=boot)
        expect_error(n_year_estimates(obj, num_years_to_estimate=nyearest,
                                      population_size = 35e5), msg)
    }

    expect_err(prevsim, 10, '2004-01-01', 8, cure=5, boot=10, nyearest=11)
    expect_err(prevsim, 5, '2003-02-23', 8, cure=5, boot=10, nyearest=8)
    expect_err(prevsim, 5, '2005-02-23', 3, cure=5, boot=10, nyearest=6)
})