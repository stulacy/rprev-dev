#' Estimate point prevalence at an index date.
#'
#' Point prevalence at a specific index date is estimated using contributions to
#' prevalence from both available registry data, and from Monte Carlo
#' simulations of the incidence and survival process, as outlined by Crouch et
#' al (2004) (see References).
#'
#' The most important parameter is \code{num_years_to_estimate}, which governs
#' the number of previous years of data to use when estimating the prevalence at
#' the index date. If this parameter is greater than the number of years of
#' known incident cases available in the supplied registry data (specified with
#' argument \code{num_registry_years}), then the remaining
#' \code{num_years_to_estimate - num_registry_years} years of incident data will
#' be simulated using Monte Carlo simulation.
#'
#' The larger \code{num_years_to_estimate}, the more accurate the prevalence
#' estimate will be, provided an adequate survival model can be fitted to the
#' registry data. It is therefore important to provide as much clean registry
#' data as possible.
#'
#' Simulated cases are marked with age and sex to enable agreement with
#' population survival data where a cure model is used, and calculation of the
#' posterior distributions of each.
#'
#' @param form Formula where the LHS is represented by a standard \code{Surv}
#'   object, and the RHS has three special function arguments: \code{age}, the
#'   column where age is located; \code{sex}, the column where sex is located;
#'   \code{entry}, the column where dates of entry to the registry are located;
#'   and \code{event}, the column where event dates are located.
#'
#'   This formula is used in the following way:
#'
#'   \code{Surv(time, status) ~ age(age_column_name) + sex(sex_column_name) +
#'   entry(entry_column_name) + event(event_column_name)}
#'
#'   Using the supplied \code{prevsim} dataset, it is therefore called with:
#'
#'   \code{Surv(time, status) ~ age(age) + sex(sex) + entry(entrydate) +
#'   event(eventdate)}
#' @param data A data frame with the corresponding column names provided in
#'   \code{form}.
#' @param num_years_to_estimate Number of years of data to consider when
#'   estimating point prevalence; multiple values can be specified in a vector.
#'   If any values are greater than the number of years of registry data
#'   available before \code{index_date}, incident cases
#'   for the difference will be simulated.
#' @param index_date The date at which to estimate point prevalence. Defaults to the
#' latest registry entry date.
#' @param num_reg_years The number of years of the registry for which incidence
#'   is to be calculated. Defaults to using all available complete years. Note
#'   that if more registry years are supplied than the number of years to
#'   estimate prevalence for, the survival data from the surplus registry years
#'   are still involved in the survival model fitting.
#' @param cure Integer defining cure model assumption for the calculation (in
#'   years). A patient who has survived beyond the cure time has a probability
#'   of surviving derived from the mortality rate of the general population.
#' @param N_boot Number of bootstrapped calculations to perform.
#' @param population_size Integer corresponding to the size of the population at
#'   risk.
#' @param proportion The population ratio to estimate prevalence for.
#' @param level Double representing the desired confidence interval width.
#' @param population_data A dataframe that must contain the columns \code{age},
#'   \code{rate}, and \code{sex}, where each row is the mortality rate for a
#'   person of that age and sex. Ideally, age ranges from [0, 100]. Defaults to
#'   the supplied data; see \code{\link{UKmortality}} for the format required
#'   for custom datasets.
#' @param precision Integer representing the number of decimal places required.
#' @param n_cores Number of CPU cores to run the fitting of the bootstrapped
#'   survival models. Defaults to 1; multi-core functionality is provided by the
#'   \code{doParallel} package.
#' @param start \strong{Deprecated: Use \code{index_date} instead and specify the
#' number of years of registry data to use with \code{num_reg_years}.}
#' Date from which incident cases are included in the format YYYY-MM-DD. Defaults
#' to the earliest entry date. This value is now inferred by counting back
#' \code{num_reg_years} years of registry data from the \code{index_date}. and
#' @return An S3 object of class \code{prevalence} with the following
#'   attributes: \item{estimates}{Estimated prevalence at the index date for
#'   each of the years in \code{num_years_to_estimate}.} \item{simulated}{A list
#'   containing items related to the simulation of prevalence contributions, see
#'   \code{\link{prevalence_simulated}}}. \item{counted}{Contributions to
#'   prevalence from each of the supplied registry years, see
#'   \code{\link{prevalence_counted}}.} \item{start_date}{The starting date of
#'   the registry data included in the estimation.} \item{index_date}{The index
#'   date at which the point prevalence was calculated for.}
#'   \item{known_inc_rate}{The known incidence rate for years included in the
#'   registry.} \item{nregyears}{Number of years of registry data that were
#'   used.} \item{nbootstraps}{The number of bootstrapped survival models fitted
#'   during the calculation.} \item{pval}{The p-value resulting from the
#'   chi-square test between the simulated and counted prevalent cases for the
#'   years of registry data available.} \item{y}{The Surv object used as the
#'   response in the survival modeling.} \item{means}{The covariate means from
#'   the data.}
#'
#' @references Crouch, Simon, et al. "Determining disease prevalence from
#'   incidence and survival using simulation techniques." Cancer epidemiology
#'   38.2 (2014): 193-199.
#' @examples
#' data(prevsim)
#'
#' \dontrun{
#' prevalence(Surv(time, status) ~ age(age) + sex(sex) + entry(entrydate) + event(eventdate),
#'            data=prevsim, num_years_to_estimate = c(5, 10), population_size=1e6,
#'            index_date = '2013-09-01', num_reg_years = 8,
#'            cure = 5)
#'
#' prevalence(Surv(time, status) ~ age(age) + sex(sex) + entry(entrydate) + event(eventdate),
#'            data=prevsim, num_years_to_estimate = 5, population_size=1e6)
#'
#' # Run on multiple cores
#' prevalence(Surv(time, status) ~ age(age) + sex(sex) + entry(entrydate) + event(eventdate),
#'            data=prevsim, num_years_to_estimate = c(3,5,7), population_size=1e6, n_cores=4)
#' }
#'
#' @export
#' @family prevalence functions
prevalence <- function(form, data, num_years_to_estimate, index_date=NULL,
                       num_reg_years=NULL, cure=10,
                       N_boot=1000,
                       population_size=NULL, proportion=100e3,
                       level=0.95, population_data=NULL,
                       precision=2, n_cores=1,
                       start=NULL) {

    # Extract required column names from formula
    spec <- c('age', 'sex', 'entry', 'event')
    terms <- terms(form, spec)
    special_indices <- attr(terms, 'specials')

    if (any(sapply(special_indices, is.null)))
        stop("Error: provide function terms for age, sex, entry date, and event date.")

    v <- as.list(attr(terms, 'variables'))[-1]
    var_names <- unlist(lapply(special_indices, function(i) v[i]))

    age_var <- .extract_var_name(var_names$age)
    sex_var <- .extract_var_name(var_names$sex)
    entry_var <- .extract_var_name(var_names$entry)
    event_var <- .extract_var_name(var_names$event)

    # Extract survival formula
    response_index <- attr(terms, 'response')
    resp <- v[response_index][[1]]
    survobj <- with(data, eval(resp))

    # Other covariates
    non_covariate_inds <- c(response_index, unlist(special_indices))
    covar_names <- as.list(attr(terms, 'variables'))[-1][-non_covariate_inds]  # First -1 to remove 'list' entry
    if (length(covar_names) > 0)
        stop("Error: functionality isn't currently provided for additional covariates.")

    # Remove when deprecate
    if ((!is.null(start) & is.null(index_date))) {
        warning("'start' parameter is deprecated and will be removed in future releases. Specify index_date and num_reg_years instead.")
        start_date <- start
        num_reg_years_new <- if (is.null(num_reg_years)) floor(as.numeric(difftime(max(data[, entry_var]), start_date) / DAYS_IN_YEAR)) else num_reg_years
        foo <- determine_yearly_endpoints(start_date, num_reg_years_new)
        index_date = foo[length(foo)]
    } else {
        if (is.null(index_date)) {
            index_date <- max(data[, entry_var])
            message("Index date not provided, using last entry date instead: ", index_date)
        }

        num_reg_years_new <- if (is.null(num_reg_years)) floor(as.numeric(difftime(index_date, min(data[, entry_var])) / DAYS_IN_YEAR)) else num_reg_years
        start_date <- as.character(as.Date(index_date) - num_reg_years_new * DAYS_IN_YEAR)
    }

    # Calculate observed prevalence for 1:num_registry_years
    # Remove calls to num_registry_years and start when update
    prev_count <- prevalence_counted(data[, entry_var],
                                     data[, event_var],
                                     survobj[, 2],
                                     index_date=index_date,
                                     num_reg_years=num_reg_years_new)


    if (num_reg_years_new > max(num_years_to_estimate)) {
        msg <- paste("More years of registry data provided than prevalence is to be estimated from.",
                     "Prevalence will be counted from the available registry data with no simulation required.")
        message(msg)
        prev_sim <- NA
        inc_rate <- rev(yearly_incidence(data[, entry_var], start_date=start_date, num_years=num_reg_years_new))
    } else {
        # Calculate simulated prevalence for 1:max(num_years_to_estimate)
        prev_sim <- prevalence_simulated(survobj, data[, age_var], data[, sex_var], data[, entry_var],
                                         max(num_years_to_estimate), index_date,
                                         num_reg_years_new, cure=cure, N_boot=N_boot,
                                         population_data=population_data, n_cores=n_cores)
        inc_rate <- prev_sim$known_inc_rate
        prev_sim$known_inc_rate <- NULL
    }


    # Calculate CIs and iterate for every year of interest
    names <- sapply(num_years_to_estimate, function(x) paste('y', x, sep=''))
    estimates <- lapply(setNames(num_years_to_estimate, names), .point_estimate,
                        prev_sim, prev_count, num_reg_years_new, population_size, level=level, precision=precision,
                        proportion=proportion)

    object <- list(estimates=estimates, simulated=prev_sim,
                   counted=prev_count, start_date=start_date,
                   index_date=index_date, known_inc_rate=inc_rate,
                   nregyears=num_reg_years_new, proportion=proportion)

    # Calculate covariate means and save
    mean_df <- data[, c(age_var, sex_var)]
    mean_df <- apply(mean_df, 2, as.numeric)
    object$means <- colMeans(mean_df)
    object$y <- survobj

    object$pval <- test_prevalence_fit(object)

    attr(object, 'class') <- 'prevalence'
    object
}

#' Count prevalence from registry data.
#'
#' Counts contribution to prevalence at a specific index from each year of a
#' registry. A person is included as contributing to disease prevalence if they
#' are incident within the specified time-span, and are either alive or censored
#' at the index date. The rationale for including censored cases in prevalence
#' estimation is that such cases have typically been lost to follow-up, and are
#' often more likely to have been alive at the index date than not.
#'
#' @inheritParams prevalence
#' @param entry Vector of diagnosis dates for each patient in the registry in
#'   the format YYYY-MM-DD.
#' @param eventdate Vector of dates corresponding to the indicator variable in
#'   the format YYYY-MM-DD.
#' @param status Vector of binary values indicating if an event has occurred for
#'   each patient in the registry. \code{entry}, \code{eventdate}, and
#'   \code{status} must all have the same length.
#' @return A vector of length equal to the number of complete years of registry data
#'   before \code{index_date} (given by \code{num_reg_years}), representing the number of
#'   incident cases in the corresponding year that contribute to the prevalence
#'   at the index date.
#' @examples
#' data(prevsim)
#'
#' prevalence_counted(prevsim$entrydate,
#'                    prevsim$eventdate,
#'                    prevsim$status)
#'
#' prevalence_counted(prevsim$entrydate,
#'                    prevsim$eventdate,
#'                    prevsim$status,
#'                    index_date="2012-01-30")
#'
#' @export
#' @family prevalence functions
prevalence_counted <- function(entry, eventdate, status, index_date=NULL, num_reg_years=NULL, start=NULL ) {

    entry <- as.Date(entry)

    if (length(unique(c(length(entry), length(eventdate), length(status)))) > 1)
        stop("Error: entry, eventdate, and status must all have the same length.")

    # Remove when deprecate
    if (!is.null(start) & is.null(index_date)) {
        warning("'start' parameter is deprecated and will be removed in future releases. Specify index_date and num_reg_years instead.")
        start_date <- if (is.null(start)) min(entry) else start
        num_reg_years_new <- if (is.null(num_reg_years)) floor(as.numeric(difftime(max(entry), start_date) / DAYS_IN_YEAR)) else num_reg_years
        foo <- determine_yearly_endpoints(start_date, num_reg_years_new)
        index_date = foo[length(foo)]
    } else {

        if (is.null(index_date)) {
            index_date <- max(entry)
            message("Index date not provided, using last entry date instead: ", index_date)
        }

        num_reg_years_new <- if (is.null(num_reg_years)) floor(as.numeric(difftime(index_date, min(entry)) / DAYS_IN_YEAR)) else num_reg_years
        start_date <- as.character(as.Date(index_date) - num_reg_years_new * DAYS_IN_YEAR)
    }


    registry_years <- determine_yearly_endpoints(start_date, num_reg_years_new)
    # Need no NAs for this!
    clean <- complete.cases(entry) & complete.cases(eventdate) & complete.cases(status)
    entry <- entry[clean]
    eventdate <- eventdate[clean]
    status <- status[clean]

    status_at_index <- ifelse(eventdate > index_date, 0, status)

    per_year <- yearly_incidence(entry, start_date=start_date, num_years=num_reg_years_new)
    num_cens <- vapply(seq(num_reg_years_new), function(i)
                            sum(status_at_index[entry >= registry_years[i] & entry < registry_years[i + 1]]),
                       numeric(1))
    per_year - num_cens
}

#' Estimate prevalence using Monte Carlo simulation.
#'
#' Estimates prevalent cases at a specific index date by use of Monte Carlo
#' simulation. Simulated cases are marked with age and sex to enable agreement
#' with population survival data where a cure model is used, and calculation of
#' the posterior distributions of each.
#'
#' @inheritParams prevalence
#' @param survobj \code{Surv} object from \code{survival} package. Currently
#'   only right censoring is supported.
#' @param age A vector of ages from the registry.
#' @param sex A vector of sex, encoded as 0 and 1 for males and females
#'   respectively.
#' @param entry A vector of entry dates into the registry, in the format
#'   YYYY-MM-DD.
#' @return A list with the following attributes:
#'   \item{mean_yearly_contributions}{A vector of length
#'   \code{num_years_to_estimate}, representing the average number of prevalent
#'   cases subdivided by year of diagnosis across each bootstrap iteration.}
#'   \item{posterior_age}{Posterior distributions of age, sampled at every
#'   bootstrap iteration.} \item{yearly_contributions}{Total simulated prevalent
#'   cases from every bootstrapped sample.} \item{pop_mortality}{Population
#'   survival rates in the format of a list, stratified by sex.}
#'   \item{nbootstraps}{Number of bootstrapped samples used in the prevalence
#'   estimation.} \item{coefs}{The bootstrapped Weibull coefficients used by the
#'   survival models.} \item{full_coefs}{The Weibull coefficients from a model
#'   fitted to the full dataset.}
#' @examples
#' data(prevsim)
#'
#' \dontrun{
#' prevalence_simulated(Surv(prevsim$time, prevsim$status), prevsim$age,
#'                      prevsim$sex, prevsim$entrydate,
#'                      num_years_to_estimate = 10,
#'                      index_date = "2013-09-01",
#'                      num_reg_years = 8, cure = 5)
#'
#' prevalence_simulated(Surv(prevsim$time, prevsim$status), prevsim$age,
#'                      prevsim$sex, prevsim$entrydate,
#'                      num_years_to_estimate = 5,
#'                      index_date="2009-01-01",
#'                      num_reg_years=5)
#'
#' # The program can be run using parallel processing.
#' prevalence_simulated(Surv(prevsim$time, prevsim$status), prevsim$age,
#'                      prevsim$sex, prevsim$entrydate,
#'                      num_years_to_estimate = 10,
#'                      index_date="2013-01-01",
#'                      num_reg_years=8, n_cores=4)
#' }
#'
#' @importFrom utils data
#' @import stats
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @export
#' @family prevalence functions
prevalence_simulated <- function(survobj, age, sex, entry, num_years_to_estimate,
                                 index_date, num_reg_years, cure=10, start=NULL,
                                 N_boot=1000,
                                 population_data=NULL, n_cores=1) {

    cure_days <- cure * DAYS_IN_YEAR

    sex <- as.factor(sex)
    if (length(levels(sex)) > 2)
        stop("Error: function can't currently function with more than two levels of sex.")

    # Calculate population survival rates for each sex in dataset
    if (is.null(population_data)) {
        utils::data('UKmortality', envir=environment())
        #assign('population_data', UKmortality)
        population_data <- get('UKmortality', envir=environment())

        #population_data <- UKmortality
    } else {
        # Obtain population data, and ensure it has the correct columns
        req_pop_names <- c('rate', 'age', 'sex')
        if (!all(sapply(req_pop_names, function(x) x %in% names(population_data)))) {
           stop("Error: the supplied population data frame must contain columns 'rate', 'age', 'sex'.")
        }
    }

    if (!is.null(start)) {
        warning("'start' parameter is deprecated and will be removed in future versions. Please use specify an index date using 'index_date' instead.")
    }

    start <- as.character(as.Date(index_date) - num_reg_years * DAYS_IN_YEAR)

    population_data$sex <- as.factor(population_data$sex)
    if (!all(levels(sex) %in% levels(population_data$sex)))
        stop("Error: the same levels must be present in both the population and registry data. '0' and '1' by default where male is 0.")

    surv_functions <- lapply(setNames(levels(sex), levels(sex)),
                             function(x) population_survival_rate(rate ~ age, data=subset(population_data, sex==x)))

    # Calculate bootstrapped weibull coefficients for registry data
    df <- data.frame(time=survobj[, 1], status=survobj[, 2], age=age, sex=sex, entry=entry)
    df <- df[as.character(df$entry) >= start, ]

    # Specify whether to include sex as a survival variable or not, this should be included within the formula!
    req_covariate <- if (length(levels(sex)) == 1) 'age' else 'age + sex'
    surv_form <- as.formula(paste('Surv(time, status)', '~',
                                  req_covariate))
    wb_boot <- .registry_survival_bootstrapped(surv_form, df, N_boot, n_cores=n_cores)
    wb_boot <- wb_boot[sample(nrow(wb_boot)), ]

    # Calculate maximum yearly incidence here as the max of the maximum for each subgroup with an additional bit of leeway
    subgroup_max_incidence <- sapply(levels(df$sex), function(x) {
        sub_data <- df[df$sex==x, ]
        max(yearly_incidence(as.character(sub_data$entry), start_date=start, num_years=num_reg_years))
    })

    # Run the prevalence estimator for each subgroup
    results <- lapply(levels(df$sex), function(x) {
        sub_data <- df[df$sex==x, ]
        .prevalence_subgroup(sub_data$age, as.character(sub_data$entry), start, wb_boot, num_reg_years,
                             surv_functions[[x]], cure_days, as.numeric(x), num_years_to_estimate,
                             include_sex=length(levels(df$sex)) == 2)
    })

    # Combine results if have more than 1 sex subgroup
    if (length(results) > 1) {
        # This is ugly but will work provided there aren't more than 2 sexs specified for, which is guarded
        # against anyway
        by_year_samples <- results[[1]]$cases + results[[2]]$cases
        post_age_dist <- lapply(results, function(x) x$post)
        fix_rate <- rev(results[[1]]$fix + results[[2]]$fix)
    } else {
        by_year_samples <- results[[1]]$cases
        post_age_dist <- results[[1]]$post
        fix_rate <- results[[1]]$fix
    }
    by_year_avg <- rowMeans(by_year_samples)

    # Fit weibull model to full data
    full_data_trans <- .transform_registry_data(surv_form, df)
    full_coefs <- .fit_weibull(full_data_trans)
    full_coefs[length(full_coefs)] <- exp(full_coefs[length(full_coefs)])  # As this .fit_weibull skips the reverse log transform

    prev_out <- list(mean_yearly_contributions=by_year_avg, posterior_age=post_age_dist,
                     yearly_contributions=by_year_samples, known_inc_rate=fix_rate,
                     pop_mortality=surv_functions, nbootstraps=N_boot, coefs=wb_boot, full_coefs=full_coefs)
    prev_out
}


.prevalence_subgroup <- function(prior_age_d, entry, start, wboot, nregyears, survfunc,
                                 cure_days, sex, nprevyears, include_sex) {
    fix_rate_rev <- rev(yearly_incidence(entry, start_date=start, num_years=nregyears))
    mean_rate <- mean(fix_rate_rev)

    #  This is the new implementation of calculating the yearly predicted prevalence
    yearly_rates = lapply(1:nprevyears, .yearly_prevalence, wboot, mean_rate, nregyears, fix_rate_rev,
                          prior_age_d, survfunc, cure_days, sex, include_sex)
    # Unflatten by_year samples
    by_year_samples = do.call(rbind, lapply(yearly_rates, function(x) x$cases))
    post_age_dist = lapply(yearly_rates, function(x) x$post)
    return(list(cases=by_year_samples, post=post_age_dist, fix=fix_rate_rev))
}


.yearly_prevalence <- function(year, bootwb, meanrate, nregyears, fixrate, prior, dailysurv, cure_days, sex, include_sex) {
    # Run the bootstrapping to obtain posterior distributions and # cases for this year
    post_results = apply(bootwb, 1, .post_age_bs, meanrate, nregyears, year-1, fixrate[year], prior, dailysurv,
                         cure_days, sex, inreg=year<=nregyears, include_sex=include_sex)

    # Post_age_bs returns a list with 'cases' and 'post' values for the number of cases and posterior age distribution
    # Need to flatten this into single array for boot_out and 2D array for post_age_dist
    bs_cases <- vapply(post_results, function(x) x$cases, integer(1))
    bs_post <- lapply(post_results, function(x) x$post)
    return(list(cases=bs_cases, post=bs_post))
}


.post_age_bs <- function(coefs_bs, meanrate, nregyears, year, fixrate, prior, daily_surv, cure_days,
                         sex, inreg=TRUE, include_sex=TRUE) {
    if (inreg){
        rate <- fixrate
    } else{
        rate <- max(0, rnorm(1, meanrate, sqrt(meanrate) / nregyears))
    }

    num_diag <- rpois(1, rate)

    if (num_diag > 0) {
        boot_age_dist <- sample(prior, num_diag, replace=TRUE)
        time_since_diag <- year * DAYS_IN_YEAR + runif(num_diag, 0, DAYS_IN_YEAR)

        # Combine data into a matrix
        bootstrapped_data <- matrix(c(rep(1, num_diag), boot_age_dist), nrow=num_diag)
        if (include_sex) {
            bootstrapped_data <- cbind(bootstrapped_data, rep(sex, num_diag))
        }

        is_dead <- as.logical(rbinom(num_diag, 1,
                                      1 - .prob_alive(time_since_diag, bootstrapped_data,
                                                      cure_days, boot_coefs=coefs_bs,
                                                      pop_surv_rate=daily_surv)))
        num_alive <- sum(!is_dead)
    } else {
        num_alive <- as.integer(0)
    }

    if (num_alive > 0) {
        post_age_dist <- time_since_diag[!is_dead]/DAYS_IN_YEAR + boot_age_dist[!is_dead]
    } else {
        post_age_dist <- numeric()
    }

    return(list(cases=num_alive, post=post_age_dist))
}


#' @export
print.prevalence <- function(x, ...) {
    cat(paste0("Estimated prevalence at ", x$index_date, ":\n"))
    lapply(names(x$estimates), function(item) {
        year <- strsplit(item, 'y')[[1]][2]
        abs_prev_est <- x$estimates[[item]][1]
        if (length(x$estimates[[item]]) > 1) {
            rel_prev <- x$estimates[[item]][2]
            rel_prev_est <- paste0("(", rel_prev, " per ", x$proportion, ")")
        } else {
            rel_prev_est <- NULL
        }
        cat(paste(year, "years:", abs_prev_est, rel_prev_est, '\n'))
    })
}

#' @export
summary.prevalence <- function(object, ...) {
    cat("Prevalence \n~~~~~~~~~~\n")
    print(object)
    cat("\n")

    cat("Registry Data\n~~~~~~~~~~~~~\n")
    cat("Index date:", object$index_date, "\n")
    cat("Start date:", object$start_date, "\n")
    cat("Number of years:", length(object$known_inc_rate), "\n")
    cat("Known incidence rate:\n")
    cat(object$known_inc_rate, "\n")
    cat("Counted prevalent cases:\n")
    cat(object$counted)

    if (!all(is.na(object$simulated))) {
        cat("\n\nBootstrapping\n~~~~~~~~~~~~~\n")
        cat("Iterations:", object$simulated$nbootstraps, "\n")
        cat("Posterior age distribution summary:\n")
        print(summary(unlist(object$simulated$posterior_age)))
        cat("Average simulated prevalent cases per year:\n")
        cat(round(rev(object$simulated$mean_yearly_contributions)), "\n")
        cat("P-value from chi-square test:", object$pval)
    }
}


.point_estimate <- function(years, sim, obs, num_reg_years, population_size=NULL, proportion=100e3,
                            level=0.95, precision=2) {

    # Replace simulated prevalence contributions from observed years with counted prevalence
    if (all(is.na(sim))) {
        prev_estimate <- rev(obs)
    } else {
        prev_estimate <- sim$mean_yearly_contributions
        prev_estimate[1:num_reg_years] <- rev(obs)
    }

    # Replace simulated values for observed for the years we have simulated data
    #sim$mean_yearly_contributions[1:num_reg_years] <- rev(obs)

    the_estimate <- sum(prev_estimate[1:years])
    result <- list(absolute.prevalence=the_estimate)

    if (!is.null(population_size)) {
        raw_proportion <- the_estimate / population_size
        the_proportion <- proportion * raw_proportion

        if (years <= num_reg_years) {
            se <- (raw_proportion * (1 - raw_proportion)) / population_size
        } else {
            the_samples <- sim$yearly_contributions[(num_reg_years + 1):years, , drop=FALSE]
            by_sample_estimate <- colSums(the_samples)
            the_estimate_n <- sum(prev_estimate[1:num_reg_years])
            raw_proportion_n <- the_estimate_n / population_size
            std_err_1 <- sqrt((raw_proportion_n * (1 - raw_proportion_n)) / population_size)
            std_err_2 <- sd(by_sample_estimate)/population_size
            se <- std_err_1^2 + std_err_2^2
        }

        z_level <- qnorm((1+level)/2)
        CI <- z_level * sqrt(se) * proportion

        # Setup labels for proportion list outputs
        proportion_unit <- if (proportion / 1e6 >= 1) 'M' else {
                                  if (proportion / 1e3 >= 1) 'K' else ''
                              }
        proportion_val <- if (proportion / 1e6 >= 1) proportion / 1e6 else {
                                  if (proportion / 1e3 >= 1) proportion / 1e3 else proportion
                            }
        est_lab <- paste('per', proportion_val, proportion_unit, sep='')
        upper_lab <- paste(est_lab, '.upper', sep='')
        lower_lab <- paste(est_lab, '.lower', sep='')
        result[[est_lab]] <- the_proportion
        result[[upper_lab]] <- the_proportion + CI
        result[[lower_lab]] <- the_proportion - CI
    }

    lapply(result, round, precision)
}