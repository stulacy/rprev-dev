library(rprev)
context('Survival')
data(UKmortality)
data(prevsim)

test_that("population_survival_rate returns same values as before", {
    i <- 1
    expect_ref <- function(data, age=100) {
        fn <- paste("cache/survival/population_survival_rate_", i, ".rds", sep='')
        expect_equal_to_reference(population_survival_rate(rate ~ age, data, max_age=age), file=fn)
        i <<- i + 1
    }

    expect_ref(UKmortality)
    expect_ref(subset(UKmortality, sex==0))
    expect_ref(subset(UKmortality, sex==1))
    expect_ref(UKmortality, age=50)
    expect_ref(subset(UKmortality, sex==0), age=50)
    expect_ref(subset(UKmortality, sex==1), age=50)
})

test_that("population_survival_rate returns doubles", {
    expect_double <- function(data, age=100) {
        expect_match(typeof(population_survival_rate(rate ~ age, data, max_age=age)), 'double')
    }

    expect_double(UKmortality)
    expect_double(subset(UKmortality, sex==0))
    expect_double(subset(UKmortality, sex==1))
    expect_double(UKmortality, age=50)
    expect_double(subset(UKmortality, sex==0), age=50)
    expect_double(subset(UKmortality, sex==1), age=50)
})

test_that("population_survival_rate returns the correct number of values", {
    expect_length <- function(data, age=100) {
        expect_equal(length(population_survival_rate(rate ~ age, data, max_age=age)), age*365)
    }

    expect_length(UKmortality)
    expect_length(subset(UKmortality, sex==0))
    expect_length(subset(UKmortality, sex==1))
    expect_length(UKmortality, age=50)
    expect_length(subset(UKmortality, sex==0), age=50)
    expect_length(subset(UKmortality, sex==1), age=50)
})

test_that(".registry_survival_bootstrapped returns the same value as before", {
    set.seed(3)
    i <- 1
    expect_ref <- function(data, straps) {
        form = Surv(time, status) ~ age + sex
        fn <- paste("cache/survival/registry_survival_bootstrapped_", i, ".rds", sep='')
        expect_equal_to_reference(.registry_survival_bootstrapped(form, data, straps), file=fn)
        i <<- i + 1
    }

    expect_ref(prevsim, 10)
    expect_ref(prevsim, 50)
    expect_ref(prevsim, 100)

    expect_ref(subset(prevsim, sex==0), 10)
    expect_ref(subset(prevsim, sex==0), 50)
    expect_ref(subset(prevsim, sex==0), 100)

    expect_ref(subset(prevsim, sex==1), 10)
    expect_ref(subset(prevsim, sex==1), 50)
    expect_ref(subset(prevsim, sex==1), 100)
})

test_that(".registry_survival_bootstrapped with a thousand bootstraps returns the same value as before", {
    skip_on_cran()
    skip("too slow")

    set.seed(3)
    i <- 1
    expect_ref <- function(data, straps) {
        form <- Surv(time, status) ~ age + sex
        fn <- paste("cache/survival/registry_survival_bootstrapped_1000_", i, ".rds", sep='')
        expect_equal_to_reference(.registry_survival_bootstrapped(form, data, straps), file=fn)
        i <<- i + 1
    }

    expect_ref(prevsim, 1000)
    expect_ref(subset(prevsim, sex==0), 1000)
    expect_ref(subset(prevsim, sex==1), 1000)
})

test_that(".transform_registry_data returns same values as before", {
    i <- 1
    expect_ref <- function(form, data) {
        fn <- paste("cache/survival/transform_registry_data_", i, ".rds", sep='')
        expect_equal_to_reference(.transform_registry_data(form, data), fn)
        i <<- i + 1
    }

    expect_ref(Surv(time, status) ~ age, prevsim)
    expect_ref(Surv(time, status) ~ age + sex, prevsim)
    expect_ref(Surv(time, status) ~ age + sex, subset(prevsim, sex==0))
    expect_ref(Surv(time, status) ~ sex, subset(prevsim, sex==1))
})

test_that(".transform_registry_data returns double", {
    expect_double <- function(form, data) {
        expect_match(typeof(.transform_registry_data(form, data)), 'double')
    }

    expect_double(Surv(time, status) ~ age, prevsim)
    expect_double(Surv(time, status) ~ age + sex, prevsim)
    expect_double(Surv(time, status) ~ age + sex, subset(prevsim, sex==0))
    expect_double(Surv(time, status) ~ sex, subset(prevsim, sex==1))
})

test_that(".transform_registry_data returns correct dimensions", {
    expect_dim <- function(form, data, dim) {
        expect_equal(dim(.transform_registry_data(form, data)), dim)
    }

    expect_dim(Surv(time, status) ~ 1, prevsim, c(nrow(prevsim), 3))
    expect_dim(Surv(time, status) ~ age, prevsim, c(nrow(prevsim), 4))
    expect_dim(Surv(time, status) ~ age + sex, prevsim, c(nrow(prevsim), 5))
    expect_dim(Surv(time, status) ~ age + sex, subset(prevsim, sex==0), c(nrow(prevsim[prevsim$sex==0, ]), 5))
    expect_dim(Surv(time, status) ~ sex, subset(prevsim, sex==1), c(nrow(prevsim[prevsim$sex==1, ]), 4))
})

test_that(".calculate_coefficients returns the same values as before", {
    set.seed(17)
    i <- 1
    expect_ref <- function(form, data, nobs, ncoef) {
        fn <- paste("cache/survival/calculate_coefficients_", i, ".rds", sep='')
        trans_data <- .transform_registry_data(form, data)
        expect_equal_to_reference(.calculate_coefficients(trans_data, nobs, ncoef), fn)
        i <<- i + 1
    }

    expect_ref(Surv(time, status) ~ age, prevsim, nrow(prevsim), 3)
    expect_ref(Surv(time, status) ~ age+sex, prevsim, nrow(prevsim), 4)
    expect_ref(Surv(time, status) ~ age, prevsim[1:500, ], 500, 3)
    expect_ref(Surv(time, status) ~ age+sex, prevsim[1:500, ], 500, 4)
})

test_that(".calculate_coefficients returns NAs when applicable", {
    set.seed(17)
    expect_NA <- function(form, data, n_obs, n_coef) {
        trans_data <- .transform_registry_data(form, data)
        expect_equal(.calculate_coefficients(trans_data, n_obs, n_coef), rep(NA, n_coef))
    }

    negative_times <- prevsim
    negative_times[sample(nrow(prevsim), 500), 'time'] <- rnorm(500, -50, 10)
})


test_that(".calculate_bootstrapped_coefficients returns the same values as before", {
    set.seed(17)
    i <- 1
    expect_ref <- function(form, data, N) {
        fn <- paste("cache/survival/calculate_bootstrapped_coefficients_", i, ".rds", sep='')
        trans_data <- .transform_registry_data(form, data)
        expect_equal_to_reference(.calculate_bootstrapped_coefficients(trans_data, N), fn)
        i <<- i + 1
    }

    expect_ref(Surv(time, status) ~ age, prevsim, 10)
    expect_ref(Surv(time, status) ~ age + age, prevsim, 10)
    expect_ref(Surv(time, status) ~ sex, prevsim, 10)
    expect_ref(Surv(time, status) ~ age, prevsim, 50)
    expect_ref(Surv(time, status) ~ age + age, prevsim, 50)
    expect_ref(Surv(time, status) ~ sex, prevsim, 50)
})

test_that(".fit_weibull returns the same values as before", {
    set.seed(17)
    i <- 1
    expect_ref <- function(form, data) {
        fn <- paste("cache/survival/fit_weibull_", i, ".rds", sep='')
        trans_data <- .transform_registry_data(form, data)
        expect_equal_to_reference(.fit_weibull(trans_data), fn)
        i <<- i + 1
    }

    expect_ref(Surv(time, status) ~ age, prevsim)
    expect_ref(Surv(time, status) ~ age + age, prevsim)
    expect_ref(Surv(time, status) ~ sex, prevsim)
    expect_ref(Surv(time, status) ~ age, prevsim[1:500, ])
    expect_ref(Surv(time, status) ~ age + age, prevsim[1:500, ])
    expect_ref(Surv(time, status) ~ sex, prevsim[1:500, ])
})

test_that(".row_any_error correctly returns NA for a variety of errors", {
    expect_rows <- function(input, output) {
        expect_equal(.row_any_error(input), output)
    }

    expect_rows(matrix(seq(25), ncol=5), rep(F, 5))

    matrix <- matrix(seq(16), nrow=4)

    matrix1 <- matrix
    matrix1[1, 2] <- NA
    expect_rows(matrix1, c(T, F, F, F))

    matrix2 <- matrix
    matrix2[3, 1] <- NaN
    expect_rows(matrix2, c(F, F, T, F))

    matrix3 <- matrix
    matrix3[3, 1] <- NA
    matrix3[4, 2] <- NA
    expect_rows(matrix3, c(F, F, T, T))

    matrix4 <- matrix
    matrix4[1, 2] <- NaN
    matrix4[3, 4] <- 542
    expect_rows(matrix4, c(T, F, T, F))

    matrix5 <- matrix
    matrix5[1, 2] <- NA
    matrix5[2, 1] <- NaN
    matrix5[3, 4] <- 542
    matrix5[4, 2] <- 223
    expect_rows(matrix5, c(T, T, T, F))

    matrix6 <- matrix
    matrix6[1, 2] <- NaN
    matrix6[2, 1] <- 827082
    matrix6[3, 4] <- NA
    matrix6[4, 4] <- 3e3
    expect_rows(matrix6, c(T, F, T, T))
})

test_that(".prob_alive provides the correct probability of death", {
    set.seed(3)
    i <- 1
    expect_ref <- function(time, data, cure_time, coefs, pop_rate, max_age=100) {
        fn <- paste("cache/survival/prob_alive_", i, ".rds", sep='')
        expect_equal_to_reference(.prob_alive(time, data, cure_time, coefs, pop_rate, max_age), file=fn)
        i <<- i + 1
    }

    data = prevsim
    time = 5000
    cure_time = 9000
    max_age = 100

    mat_data <- matrix(c(rep(1, nrow(data)), data$age), nrow=nrow(data))
    trans_data <- .transform_registry_data(Surv(time, status) ~ age, data)
    boot_coefs = .calculate_bootstrapped_coefficients(trans_data, 5)
    pop_surv_rate = population_survival_rate(rate ~ age, UKmortality, max_age=max_age)

    # Test with different coefficients
    expect_ref(time, mat_data, cure_time, boot_coefs[3, ], pop_surv_rate, max_age)
    expect_ref(time, mat_data, cure_time, boot_coefs[4, ], pop_surv_rate, max_age)

    # Test with time > cure_time
    cure_time <- 3500
    expect_ref(time, mat_data, cure_time, boot_coefs[1, ], pop_surv_rate, max_age)

    # Test with sex included
    mat_data <- matrix(c(rep(1, nrow(data)), data$age, data$sex), nrow=nrow(data))
    trans_data <- .transform_registry_data(Surv(time, status) ~ age + sex, data)
    boot_coefs = .calculate_bootstrapped_coefficients(trans_data, 5)

    expect_ref(time, mat_data, cure_time, boot_coefs[3, ], pop_surv_rate, max_age)
    cure_time <- 6000
    expect_ref(time, mat_data, cure_time, boot_coefs[5, ], pop_surv_rate, max_age)

    # With only males
    trans_data <- .transform_registry_data(Surv(time, status) ~ age + sex, data)
    boot_coefs = .calculate_bootstrapped_coefficients(trans_data, 5)

    male_data <- prevsim[prevsim$sex==0, ]
    mat_data <- matrix(c(rep(1, nrow(male_data)), male_data$age, male_data$sex), nrow=nrow(male_data))
    pop_surv_rate = population_survival_rate(rate ~ age, UKmortality[UKmortality$sex==0, ], max_age=max_age)

    # With only females
    trans_data <- .transform_registry_data(Surv(time, status) ~ age + sex, data)
    boot_coefs = .calculate_bootstrapped_coefficients(trans_data, 5)

    female_data <- prevsim[prevsim$sex==1, ]
    mat_data <- matrix(c(rep(1, nrow(female_data)), female_data$age, female_data$sex), nrow=nrow(female_data))
    pop_surv_rate = population_survival_rate(rate ~ age, UKmortality[UKmortality$sex==1, ], max_age=max_age)

})
