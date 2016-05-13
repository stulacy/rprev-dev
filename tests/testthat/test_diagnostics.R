library(prevR)
library(survival)
library(rms)
library(abind)
context('Diagnostics')
data(prevsim)

test_that("sim_check returns same values as before", {
    expect_ref <- function(data, N_sim) {
        fn <- 'cache/diagnostics/sim_check.rds'
        set.seed(17)
        expect_equal_to_reference(sim_check(incidence(prevsim$entrydate), N_sim = 10), file=fn)
    }
    set.seed(17)
    expect_ref(incidence(prevsim$entrydate), N_sim = 10)
})

test_that("sim_check with 100000 simulations returns same values as before", {
    skip_on_cran()
    skip("too slow")
    expect_ref <- function(data) {
        fn <- 'cache/diagnostics/sim_check_100000.rds'
        set.seed(18)
        expect_equal_to_reference(sim_check(incidence(prevsim$entrydate)), file=fn)
    }
    set.seed(18)
    expect_ref(incidence(prevsim$entrydate))
})

test_that("sim_check returns doubles", {
    expect_double <- function(data, N_sim) {
        expect_match(typeof(sim_check(incidence(prevsim$entrydate), N_sim = 10)), 'double')
    }

    expect_double(incidence(prevsim$entrydate), N_sim = 10)
})

test_that("sim_check returns no NAs", {
    expect_NA <- function(data, N_sim) {
        expect_equal(any(is.na(sim_check(incidence(prevsim$entrydate), N_sim = 10))), FALSE)
    }

    expect_NA(incidence(prevsim$entrydate), N_sim = 10)
})

test_that("sim_check returns the correct number of values", {
    expect_length <- function(data, N_sim) {
        expect_equal(length(sim_check(incidence(prevsim$entrydate), N_sim = 10)),2)
    }

    expect_length(sim_check(incidence(prevsim$entrydate), N_sim = 10))
})

test_that("functional_form_age returns a list", {
    expect_list <- function(data) {
        expect_match(typeof(functional_form_age(Surv(time, status) ~ age(age), data)), 'list')
    }

    expect_list(prevsim)
})

test_that("prev_chisq returns same values as before", {
    expect_ref <- function(data) {
        fn <- 'cache/diagnostics/prev_chisq.rds'
        set.seed(17)
        prevalence_object <- prevalence(Surv(time, status) ~ sex(sex) + age(age) + entry(entrydate),
                                        data, num_years_to_estimate=10,
                                        start='2004-01-30', num_reg_years=9)
        expect_equal_to_reference(prev_chisq(prevalence_object), file=fn)
    }

    expect_ref(prevsim)
})

test_that("prev_chisq returns doubles", {
    expect_double <- function(data) {
        set.seed(17)
        prevalence_object <- prevalence(Surv(time, status) ~ sex(sex) + age(age) + entry(entrydate),
                                        data, num_years_to_estimate=10,
                                        start='2004-01-30', num_reg_years=9)
        expect_match(typeof(prev_chisq(prevalence_object)), 'double')
    }

    expect_double(prevsim)
})

test_that("prev_chisq returns no NAs", {
    expect_NA <- function(data) {
        set.seed(17)
        prevalence_object <- prevalence(Surv(time, status) ~ sex(sex) + age(age) + entry(entrydate),
                                        data, num_years_to_estimate=10,
                                        start='2004-01-30', num_reg_years=9)
        expect_equal(any(is.na(prev_chisq(prevalence_object))), FALSE)
    }

    expect_NA(prevsim)
})

test_that("cumulative_incidence returns same values as before", {
    set.seed(3)
    expect_ref <- function(data) {
        fn <- 'cache/diagnostics/cumulative_incidence.rds'
        expect_equal_to_reference(cumulative_incidence(prevsim$entrydate), file=fn)
    }
    expect_ref(prevsim$entrydate)
})


test_that("inspect_incidence returns same as before", {
    expect_ref <- function(data) {
        fn <- 'cache/diagnostics/inspect_incidence.rds'
        c_inc <- cumulative_incidence(data)
        expect_equal_to_reference(inspect_incidence(c_inc), file=fn)
    }
    expect_ref(prevsim$entrydate)
})

test_that("poisson_incidence_sim returns same as before", {
    skip_on_cran()
    skip("too slow")
    set.seed(3)
    expect_ref <- function(data) {
        fn <- 'cache/diagnostics/poisson_incidence_sim.rds'
        c_inc <- cumulative_incidence(data)
        expect_equal_to_reference(poisson_incidence_sim(c_inc), file=fn)
    }
    expect_ref(prevsim$entrydate)
})

test_that("boot_eg returns same as before", {
    #skip_on_cran()
    #skip("too slow")
    i <- 1
    set.seed(3)
    expect_ref <- function(data, sex) {
        fn <- paste('cache/diagnostics/boot_eg', i, '.rds', sep='')
        expect_equal_to_reference(boot_eg(form = Surv(time, status) ~ age(age) + sex(sex) + entry(entrydate),
                                                       data = data,
                                                       start = "2004-01-30",
                                                       age = 65,
                                                       sex = sex), file=fn)
        i <<- i + 1
    }
    expect_ref(prevsim, "Male")
    expect_ref(prevsim, "Female")
})