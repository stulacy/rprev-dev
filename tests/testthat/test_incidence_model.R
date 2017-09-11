library(rprev)
context('Incidence model')
data(prevsim)

test_that("rate is correctly estimated within 10%", {
    expect_rate <- function(startdate, rate, n_inds=1000) {
        entrydates <- as.Date(startdate) + cumsum(rexp(n_inds, rate))
        mod <- fit_exponential_incidence(entry ~ 1, data.frame(entry=entrydates))
        error <- abs((mod$rate - rate) / rate)
        expect_lte(error, 0.10)
    }
    expect_rate('2004-01-01', 9)
    expect_rate('2005-03-23', 0.5)
    expect_rate('1998-05-17', 0.003)
    expect_rate('1885-09-12', 100)
})

test_that("rate for stratified covariates is correctly estimated within 10%", {
    expect_rate <- function(startdate, rates, n_inds=1000) {
        # Generate stratified entry date
        df <- data.frame(entry=c(), group=c())
        for (i in seq_along(rates)) {
            entrydates <- as.Date(startdate) + cumsum(rexp(n_inds, rates[i]))
            df <- rbind(df, data.frame(entry=entrydates, group=LETTERS[i]))
        }
        # Shuffle df
        df <- df[sample(1:nrow(df), replace=F), ]

        # Fit model
        mod <- fit_exponential_incidence(entry ~ group, df)
        error <- abs((mod$rate$Freq - rates) / rates)
        expect_false(any(error > 0.10))
    }
    expect_rate('2004-01-01', c(0.52, 0.32))
    expect_rate('2005-03-23', c(5, 3))
    expect_rate('1998-05-17', c(0.7, 0.5, 0.4))
    expect_rate('1885-09-12', c(50, 40, 60))
})