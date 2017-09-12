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

test_that("test_prevalence_fit returns same values as before without error and isn't significant", {
    set.seed(17)
    prevalence_object <- prevalence("2013-01-01",
                                    num_years_to_estimate=10,
                                    data=prevsim,
                                    inc_formula=entrydate ~ sex,
                                    surv_formula=Surv(time, status) ~ sex + age,
                                    dist='weibull',
                                    death_column='eventdate',
                                    N_boot=20)
    fn <- 'cache/diagnostics/prev_pval.rds'
    expect_equal_to_reference(test_prevalence_fit(prevalence_object), file=fn)
    expect_gt(prevalence_object$pval, 0.05)
    expect_match(typeof(test_prevalence_fit(prevalence_object)), 'double')
    expect_equal(any(is.na(test_prevalence_fit(prevalence_object))), FALSE)
})
