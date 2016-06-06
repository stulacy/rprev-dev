library(rprev)
context('Prevalence')
data(prevsim)

test_that("prevalence returns same values as before", {
    set.seed(3)
    i <- 1
    expect_ref <- function(data, num_years_to_estimate, start, years, cure, boot) {
        fn <- paste('cache/prevalence/prevalence_', i, '.rds', sep='')
        expect_equal_to_reference(prevalence(Surv(time, status) ~ sex(sex) + age(age) + entry(entrydate) + event(eventdate),
                                             data, num_years_to_estimate, population_size=1e6,
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
        expect_equal_to_reference(prevalence(Surv(time, status) ~ sex(sex) + age(age) + entry(entrydate) + event(eventdate),
                                             data, num_years_to_estimate, population_size=1e6,
                                             start=start, num_reg_years=years,
                                             cure=cure, N_boot=boot), file=fn)
        i <<- i + 1
    }
    expect_ref(prevsim, 10, '2004-01-01', 9, cure=5, boot=1000)
    expect_ref(prevsim, 10, '2006-09-03', 5, cure=5, boot=1000)
})

test_that("Error is raised when passing a population data frame not set up correctly", {
    expect_poperror <- function(popdata, msg) {
        expect_error(prevalence(Surv(time, status) ~ sex(sex) + age(age) + entry(entrydate) + event(eventdate),
                                             prevsim, num_years_to_estimate=10, population_size=1e6, population_data=popdata,
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
                                prevsim, num_years_to_estimate=10, population_size = 1e6,
                                start='2004-01-01', num_reg_years=9,
                                cure=5, N_boot=10), msg)
    }
    missing_funcs <- "Error: provide function terms for age, sex, entry date, and event date."
    expect_formerror(Surv(time, status) ~ age + sex + entry + event, missing_funcs)
    expect_formerror(Surv(time, status) ~ age(age) + sex + entry + event, missing_funcs)
    expect_formerror(Surv(time, status) ~ age + sex(sex) + entry + event(eventdate), missing_funcs)
    expect_formerror(Surv(time, status) ~ age + sex + entry(entry), missing_funcs)
    expect_formerror(Surv(time, status) ~ age(age) + sex + entry(entry), missing_funcs)

    # All the functions are present but the columns aren't! i.e. entry is entrydate in prevsim dataset
    missing_cols <- "undefined columns selected"
    expect_formerror(Surv(time, status) ~ age(age) + sex(sex) + entry(entry) + event(eventdate), missing_cols)
    expect_formerror(Surv(time, status) ~ age(age) + sex(gender) + entry(entrydate) + event(eventdate), missing_cols)
    expect_formerror(Surv(time, status) ~ age(age) + sex(sex) + entry(entry) + event(event), missing_cols)
    expect_formerror(Surv(time, status) ~ age(age) + sex(gender) + entry(entrydate) + event(eventdate), missing_cols)
    expect_formerror(Surv(time, status) ~ age(age) + sex(sex) + entry(entrydate) + event(event), missing_cols)

})

test_that("Error is raised when levels for sex aren't the same in registry and population data", {
    expect_sexerror <- function(regdata, popdata, msg) {
        expect_error(prevalence(Surv(time, status) ~ sex(sex) + age(age) + entry(entrydate) + event(eventdate),
                                regdata, num_years_to_estimate=10, population_size=1e6,
                                population_data=popdata,
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

    reg_01 <- data.frame(time=time, status=status, age=agereg, sex=rbinom(10, 1, 0.5),
                         entrydate=entry, eventdate=as.Date(entry) + time)
    reg_mf <- data.frame(time=time, status=status, age=agereg, sex=sample(c('M', 'F'), 10, replace=T),
                         entrydate=entry, eventdate=as.Date(entry) + time)

    expect_sexerror(reg_01, pop_mf, missing_levels)
    expect_sexerror(reg_mf, pop_01, missing_levels)
})

test_that("Prevalence messages the user when the number of registry years are greater than the number of years asked to predict prevalence for", {
    expect_msg <- function(Nyears, start, regyears=NULL) {
        expect_message(prevalence(Surv(time, status) ~ sex(sex) + age(age) + entry(entrydate) + event(eventdate),
                                  prevsim, num_years_to_estimate=Nyears, population_size=1e6,
                                  start=start, num_reg_years=regyears,
                                  cure=5, N_boot=10))
    }
    expect_msg(Nyears=5, start='2003-01-01', regyears=10)
    expect_msg(Nyears=5, start='2003-01-01')  # Note there are 10 registry years in this dataset which are used if not provided
    expect_msg(Nyears=1, start='2011-01-01')  # However if we set registry to start at 2011 then only 2 years of registry data
})

test_that("prevalence_counted results in the same values as before", {
    set.seed(3)
    i <- 1
    expect_ref <- function(entry, event, status, start, years=NULL, index=NULL) {
        fn <- paste('cache/prevalence/countedprevalence_', i, '.rds', sep='')
        expect_equal_to_reference(prevalence_counted(entry, event, status, start=start,
                                                     num_reg_years=years), file=fn)
        i <<- i + 1
    }
    expect_ref(prevsim$entrydate, prevsim$eventdate, prevsim$status, '2004-03-05', 5)
    expect_ref(prevsim$entrydate, prevsim$eventdate, prevsim$status, '2004-03-05')
})
