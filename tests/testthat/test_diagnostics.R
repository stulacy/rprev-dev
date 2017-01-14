library(rprev)
context('Diagnostics')
data(prevsim)

test_that("test_incidence_fit returns same values as before", {
    expect_ref <- function(data, N_sim) {
        fn <- 'cache/diagnostics/sim_check.rds'
        set.seed(17)
        expect_equal_to_reference(test_incidence_fit(yearly_incidence(prevsim$entrydate), N_sim = 10), file=fn)
    }
    set.seed(17)
    expect_ref(incidence(prevsim$entrydate), N_sim = 10)
})

test_that("test_incidence_fit with 100000 simulations returns same values as before", {
    skip_on_cran()
    skip("too slow")
    expect_ref <- function(data) {
        fn <- 'cache/diagnostics/sim_check_100000.rds'
        set.seed(18)
        expect_equal_to_reference(test_incidence_fit(yearly_incidence(prevsim$entrydate)), file=fn)
    }
    set.seed(18)
    expect_ref(incidence(prevsim$entrydate))
})

test_that("test_incidence_fit returns doubles", {
    expect_double <- function(data, N_sim) {
        expect_match(typeof(test_incidence_fit(yearly_incidence(prevsim$entrydate), N_sim = 10)), 'double')
    }

    expect_double(incidence(prevsim$entrydate), N_sim = 10)
})

test_that("test_incidence_fit returns no NAs", {
    expect_NA <- function(data, N_sim) {
        expect_equal(any(is.na(test_incidence_fit(yearly_incidence(prevsim$entrydate), N_sim = 10))), FALSE)
    }

    expect_NA(incidence(prevsim$entrydate), N_sim = 10)
})

test_that("test_incidence_fit returns the correct number of values", {
    expect_length <- function(data, N_sim) {
        expect_equal(length(test_incidence_fit(yearly_incidence(prevsim$entrydate), N_sim = 10)),2)
    }

    expect_length(test_incidence_fit(incidence(prevsim$entrydate), N_sim = 10))
})

test_that("functional_form_age returns a list", {
    expect_list <- function(data) {
        expect_match(typeof(functional_form_age(Surv(time, status) ~ age, data, plot_fit=F)), 'list')
    }

    expect_list(prevsim)
})

test_that("test_prevalence_fit returns same values as before", {
    expect_ref <- function(data) {
        fn <- 'cache/diagnostics/prev_chisq.rds'
        set.seed(17)
        prevalence_object <- prevalence(Surv(time, status) ~ sex(sex) + age(age) + entry(entrydate) + event(eventdate),
                                        data, num_years_to_estimate=10, population_size=1e6,
                                        start='2004-01-30', num_reg_years=9)
        expect_equal_to_reference(test_prevalence_fit(prevalence_object), file=fn)
    }

    expect_ref(prevsim)
})

test_that("test_prevalence_fit returns doubles", {
    expect_double <- function(data) {
        set.seed(17)
        prevalence_object <- prevalence(Surv(time, status) ~ sex(sex) + age(age) + entry(entrydate) + event(eventdate),
                                        data, num_years_to_estimate=10, population_size=1e6,
                                        start='2004-01-30', num_reg_years=9)
        expect_match(typeof(test_prevalence_fit(prevalence_object)), 'double')
    }

    expect_double(prevsim)
})

test_that("test_prevalence_fit returns no NAs", {
    expect_NA <- function(data) {
        set.seed(17)
        prevalence_object <- prevalence(Surv(time, status) ~ sex(sex) + age(age) + entry(entrydate) + event(eventdate),
                                        data, num_years_to_estimate=10, population_size=1e6,
                                        start='2004-01-30', num_reg_years=9)
        expect_equal(any(is.na(test_prevalence_fit(prevalence_object))), FALSE)
    }

    expect_NA(prevsim)
})
