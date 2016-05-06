library(prevR)
library(abind)
context('Prevalence')
data(prevsim)


test_that("prevalence returns same values as before", {
    set.seed(3)
    i <- 1
    expect_ref <- function(data, N_years, start, years, cure, boot) {
        fn <- paste('cache/prevalence/prevalence_', i, '.rds', sep='')
        expect_equal_to_reference(prevalence(Surv(time, status) ~ sex(sex) + age(age) + entry(entrydate), 
                                             data, N_years,
                                             start=start, num_years=years,
                                             cure_time=cure*365, N_boot=boot), file=fn)
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
    expect_ref <- function(data, N_years, start, years, cure, boot) {
        fn <- paste('cache/prevalence/prevalence_thousand_', i, '.rds', sep='')
        expect_equal_to_reference(prevalence(Surv(time, status) ~ sex(sex) + age(age) + entry(entrydate), 
                                             data, N_years,
                                             start=start, num_years=years,
                                             cure_time=cure*365, N_boot=boot), file=fn)
        i <<- i + 1
    }
    expect_ref(prevsim, 10, '2004-01-01', 9, cure=5, boot=1000)
    expect_ref(prevsim, 10, '2006-09-03', 5, cure=5, boot=1000)
})

# TODO error testing, that if N_years is longer than a set amount it throws an error

# TODO test that population data must have correctly named columns

# TODO test that formula must have age, sex, and entry functions

# TODO test that levels in population data match those in registry

