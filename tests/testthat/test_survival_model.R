library(rprev)
context('Survival model')
data(prevsim)

test_that("Incorrectly specified inputs are correctly handled", {

    test_form <- function(form, dist) {
        build_survreg(form, prevsim, dist)
    }

    # Test that if provide columns that don't exist an error is thrown
    expect_error(test_form(Surv(stime, status) ~ 1, 'weibull'))
    expect_error(test_form(Surv(time, stat) ~ 1, 'weibull'))
    expect_error(test_form(Surv(time, status) ~ ageDiag, 'weibull'))

    # Test incorrectly specified distributions
    expect_error(test_form(Surv(time, status) ~ 1, 'exp'))
    expect_error(test_form(Surv(time, status) ~ age, 'wei'))
    expect_error(test_form(Surv(time, status) ~ sex, 'lnorm'))
    expect_error(test_form(Surv(time, status) ~ age + sex, 'llog'))
    expect_error(test_form(Surv(time, status) ~ age + sex, 'gamma'))
})

test_that("Models are build with correctly specified arguments", {

    test_mod <- function(form, dist) {
        expect_error(mod <- build_survreg(form, prevsim, dist), regexp = NA)
        expect_identical(class(mod), c('list', 'survregmin'))
        expect_true(all(is.finite(mod$coefs)))
    }

    # Test that if provide columns that don't exist an error is thrown
    test_mod(Surv(time, status) ~ 1, 'weibull')
    test_mod(Surv(time, status) ~ 1, 'lognormal')
    test_mod(Surv(time, status) ~ 1, 'exponential')

    test_mod(Surv(time, status) ~ age, 'weibull')
    test_mod(Surv(time, status) ~ age, 'lognormal')
    test_mod(Surv(time, status) ~ age, 'exponential')

    test_mod(Surv(time, status) ~ sex, 'weibull')
    test_mod(Surv(time, status) ~ sex, 'lognormal')
    test_mod(Surv(time, status) ~ sex, 'exponential')

    test_mod(Surv(time, status) ~ sex + age, 'weibull')
    test_mod(Surv(time, status) ~ sex + age, 'lognormal')
    test_mod(Surv(time, status) ~ sex + age, 'exponential')
})