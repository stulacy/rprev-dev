library(rprev)
context('Incidence')
data(prevsim)

test_that("raw_incidence returns same values as before", {
    i <- 1
    expect_ref <- function(start=NULL, years=NULL) {
        fn <- paste('cache/incidence/incidence_', i, '.rds', sep='')
        expect_equal_to_reference(raw_incidence(prevsim$entrydate, start=start, num_reg_years=years), file=fn)
        i <<- i + 1
    }
    expect_ref('2004-01-01', 9)
    expect_ref('2005-03-23', 5)
    expect_ref()
})

test_that("raw_incidence returns integers", {
    expect_integer <- function(start=NULL, years=NULL) {
        expect_match(typeof(raw_incidence(prevsim$entrydate, start=start, num_reg_years=years)), 'integer')
    }

    expect_integer('2004-01-01', 9)
    expect_integer('2005-03-23', 9)
    expect_integer()
})

test_that("raw_incidence returns no NAs", {
    expect_NA <- function(start=NULL, years=NULL) {
        expect_equal(any(is.na(raw_incidence(prevsim$entrydate, start=start, num_reg_years=years))), FALSE)
    }

    expect_NA('2004-01-01', 9)
    expect_NA('2005-03-23', 9)
    expect_NA()
})

test_that("mean_incidence_rate returns same values as before", {
    i <- 1
    expect_ref <- function(start=NULL, years=NULL) {
        fn <- paste('cache/incidence/meanir_', i, '.rds', sep='')
        rawinc <- raw_incidence(prevsim$entrydate, start=start, num_reg_years=years)
        expect_equal_to_reference(mean_incidence_rate(rawinc, population_size=35e5), file=fn)
        i <<- i + 1
    }

    expect_ref('2004-01-01', 9)
    expect_ref('2005-03-23', 5)
    expect_ref()
})

test_that("mean_incidence_rate returns a list", {
    expect_list <- function(start=NULL, years=NULL) {
        rawinc <- raw_incidence(prevsim$entrydate, start=start, num_reg_years=years)
        expect_match(typeof(mean_incidence_rate(rawinc, population_size=35e5)), 'list')
    }

    expect_list('2004-01-01', 9)
    expect_list('2005-03-23', 9)
    expect_list()
})

test_that("mean_incidence_rate returns no NAs", {
    expect_NA <- function(start=NULL, years=NULL) {
        rawinc <- raw_incidence(prevsim$entrydate, start=start, num_reg_years=years)
        expect_equal(any(is.na(mean_incidence_rate(rawinc, population_size=35e5))), FALSE)
    }

    expect_NA('2004-01-01', 9)
    expect_NA('2005-03-23', 9)
    expect_NA()
})
