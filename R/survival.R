#' Transform yearly mortality rates to daily.
#'
#' Calculates the daily mortality probability for a given population stratified
#' by age, based on their yearly mortality rates.
#'
#' @param form Formula where the LHS indicates the name of the mortality rate
#'   column, and the RHS is the column where age is located. This function
#'   assumes that the population data frame has already been stratified by sex,
#'   or any other categorical variable of interest.
#' @param data Data frame of population mortality stratified by sex and age. The
#'   following columns must be present: \code{sex}, \code{age}, and \code{rate}.
#' @param max_age Maximum age to calculate mortality for.
#' @return An estimate of the survival rate by age in days, with \code{max_age}
#'   * 365 values.
#' @examples
#' data(UKmortality)
#'
#' population_survival_rate(rate ~ age, UKmortality)
#' population_survival_rate(rate ~ age, subset(UKmortality, sex==0))
#'
#' @export
population_survival_rate <- function(form, data, max_age=100){
    # Could probably improve on this extraction of response and variable
    age_var <- as.character(form[[3]])
    rate_var <- as.character(form[[2]])

    rate <- vapply(seq(0, max_age-1),
                   function(x) mean(data[floor(data[,age_var]) == x, rate_var]),
                   numeric(1))

    a_rate <- c(2 * rate[1] - rate[2], rate, 2 * rate[max_age] - rate[max_age-1])
    base <- 183 + 365 * (1:max_age) # Where does 183/182 come from?
    base <- c(-182, 183, base)
    daily_rate <- approx(base, a_rate, 1:(max_age * 365))
    daily_rate <- daily_rate$y / 365
    cumprod(1 - daily_rate)
}

.registry_survival_bootstrapped <- function(form, data, N_boot = 1000, n_cores=1){
    data_trans <- .transform_registry_data(form, data) # df -> matrix, log transform stimes. Required for survreg.fit
    coefs <- .calculate_bootstrapped_coefficients(data_trans, N_boot, n_cores=n_cores)

    if (sum(.row_any_error(coefs)) > 1)
        warning("Error in coxph, possibly due to small number of events in bootstrap sample. Replacing with a new sample.")

    # Keep bootstrapping new samples to get non-NA coefficients
    while (sum(.row_any_error(coefs)) > 1) {
        is_error <- .row_any_error(coefs)
        coefs[is_error] <- .calculate_bootstrapped_coefficients(data_trans, sum(is_error))
    }

    # Need to transform scale back from log scale
    coefs[, ncol(coefs)] <- exp(coefs[, ncol(coefs)])
    coefs
}

.calculate_bootstrapped_coefficients <- function(data, N, n_cores=1) {
    n_obs <- nrow(data)
    n_coef <- ncol(data) - 1  # minus 2 for time + status, + 1 for shape

    if (n_cores > 1) {
        doParallel::registerDoParallel(n_cores)
        parobj <- foreach::foreach(i=1:N, .options.snow=list(preschedule=T), .packages=c('survival'), .combine='rbind',
                .inorder=F, .export='.calculate_coefficients')
        foreach::"%dopar%"(parobj, .calculate_coefficients(data, n_obs, n_coef))
    } else {
        t(replicate(N, .calculate_coefficients(data, n_obs, n_coef)))
    }
}

.calculate_coefficients <- function(data, nobs, ncoef) {
    # Helper function to bootstrap the data and fit a weibull
    # Split into a separate function to reduce repetition when called in parallel or serial
    # TODO Rename this function something more meaningful and less easily confused with calculate_bootstrapped_coefficients
    bstrap <- data[sample(1:nobs, nobs, replace=TRUE), ]
    tryCatch(.fit_weibull(bstrap),
             error=function(cond) rep(NA, ncoef)
    )
}

#' @importFrom survival survreg.fit
#' @importFrom survival survreg.control
#' @importFrom survival survreg.distributions
#' @importFrom survival Surv
.fit_weibull <- function(data) {
    # Helper function to fit weibull by calling the lower level survreg.fit function
    # Abstracts away this complex and ugly function call
    model <- survival::survreg.fit(data[, 3:ncol(data)],
                         data[, 1:2],
                         NULL, # weights
                         numeric(nrow(data)), # offset
                         NULL, # init
                         survival::survreg.control(), # controlvars
                         survival::survreg.distributions[['extreme']], # dist
                         0, # scale
                         1, # nstrat
                         0, # strata
                         NULL # parms
    )
    stats::coef(model)
}

.row_any_error <- function(matrix) {
    # Matrix is a N row matrix, this function calculates which rows of this matrix contain errors
    # Returns a logical vector of N length
    # 3 types of errors:
    #   - NA
    #   - NaN
    #   - Too large scale parameter, leads to problems with infinite log or obtaining cumulative prob
    apply(matrix, 1, function(x) any(sapply(x, function(y) is.na(y) | is.nan(y))) | x[length(x)] > 100)
}

.transform_registry_data <- function(form, data) {
    # Transforms the registry data into the format specified by survreg.fit,
    # i.e. as a matrix of values with the survival times log transformed.
    complete <- data[complete.cases(data), ]
    X <- model.matrix(form, complete)
    survobj <- with(complete, eval(form[[2]]))
    Y <- cbind(log(survobj[, 1]), survobj[, 2])
    cbind(Y, X)
}

.prob_alive <- function(time, data, cure_days, boot_coefs, pop_surv_rate, max_age=100){
    scale <- exp(boot_coefs[-length(boot_coefs)] %*% t(data)) + 0.000001  # Hack to ensure scale != 0
    shape <- 1 / boot_coefs[length(boot_coefs)]
    age <- data[, 2]  # hardcoded I know, but it will always be first after intercept
    ifelse(age*365 + time > (max_age * 365),
           0,
           ifelse(time < cure_days,
                  1 - pweibull(time, scale=scale, shape=shape),
                  (1 - pweibull(cure_days, scale=scale, shape=shape)) * pop_surv_rate[age*365 + time]/pop_surv_rate[age*365 + cure_days]))
}