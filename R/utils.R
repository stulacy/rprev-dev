print.cincidence <- function(object, ...) {
    object
}

summary.cincidence <- function(object, ...) {
    cat("Registry Data\n~~~~~~~~~~~~~\n")
    cat("Number of years:", length(object$raw_incidence), "\n")

    cat("\n\nIncidence\n~~~~~~~~~\n")

    cat("Known incidence by year:", object$raw_incidence, "\n")

    cat("Diagnoses (time since registry began):\n")
    print(summary(object$ordered_diagnoses))

    cat("Fitted smooth:\n")
    print(object$smooth)
}

print.prevalence <- function(object, ...) {
    cat("Estimated prevalence per", object$proportion, "at", object$index_date, "\n")
    lapply(names(object$estimates), function(x) {
        year <- strsplit(x, 'y')[[1]][2]
        prev_est <- object$estimates[[x]][2]
        cat(paste(year, "years:", prev_est, "\n"))
    })
}

summary.prevalence <- function(object, ...) {
    cat("Registry Data\n~~~~~~~~~~~~~\n")
    cat("Index date:", object$index_date, "\n")
    cat("Start year:", object$start_date, "\n")
    cat("Number of years:", length(object$known_inc_rate), "\n")
    cat("Known incidence rate:\n")
    cat(object$known_inc_rate, "\n")
    cat("Counted prevalent cases:\n")
    cat(object$counted)

    cat("\n\nBootstrapping\n~~~~~~~~~~~~~\n")
    cat("Iterations:", object$simulated$nbootstraps, "\n")
    cat("Posterior age distribution summary:\n")
    print(summary(object$simulated$posterior_age))
    cat("Average simulated prevalent cases per year:\n")
    cat(round(rev(object$simulated$mean_yearly_contributions)), "\n")
    cat("P-value from chi-square test:", object$pval)
}

determine_registry_years <- function(start, num_reg_years) {
    # Calculate registry years from this
    # NB: Ugly hack to not take leap years into account. Done so that tests don't throw an error, but strictly
    # should just use as.Date(start) + 365.25 * x to account for leap years
    sapply(0:num_reg_years, function(x) paste(as.numeric(strftime(start, '%Y'))+x,
                                          strftime(start, '%m-%d'), sep='-'))
}

# Extracts the variable names from a formula
# expr (language): Expression in the form entry(diagdate), where
# 'entry' is the known function name while the name inside brackets
# is the user's column name for that variable
.extract_var_name <- function(expr) {
    as.character(expr)[[2]]
}

#' Calculates bootstrapped survival probabilities from the Weibull models fitted to the
#' \code{prevalence} object.
#'
#' @param object A \code{prevalence} object.
#' @param newdata A list or dataframe with the co-variate values to calculate survival probabilities
#' for. Defaults to using the mean values from the training dataset.
#' @return An S3 object of class \code{survfit.prev} with the following attributes:
#' \item{time}{A vector of time points at which survival probability has been calculated.}
#' \item{surv}{A matrix of survival probabilities, where the rows represent a different bootstrapped
#' Weibull model, and the columns are each timepoint.}
#' \item{fullsurv}{A vector of survival probabilities for the model built on the full training set.}
#' @examples
#' data(prevsim)
#'
#' prev_obj <- prevalence(Surv(time, status) ~ age(age) + sex(sex) + entry(entrydate) + event(eventdate),
#'                        data=prevsim, num_years_to_estimate = c(5, 10), population_size=1e6,
#'                        start = "2005-09-01",
#'                        num_reg_years = 8, cure = 5)
#'
#' survobj <- survfit(prev_obj)
#'
#' survobj <- survfit(prev_obj, newdata=list(age=65, sex=0))
#'
#' @export survfit.prevalence
survfit.prevalence <- function(object, newdata=NULL) {
    if (is.null(newdata)) {
        use_df <- object$means
    } else {
        # Check names have same names as those in original data
        if (!all(sort(names(newdata)) == sort(names(object$means))))
            stop("Error: Please provide a list with the same column names as in the orignal dataset.")

        # Reorder new data to be in same order as original data
        use_df <- unlist(newdata)[names(object$means)]
    }

    # Add intercept
    use_df <- c(1, use_df)

    times <- seq(max(object$y[, 1]))
    survprobs <- apply(object$simulated$coefs, 1, .calculate_survival_probability, use_df, times)

    full_survprobs <- .calculate_survival_probability(object$simulated$full_coefs, use_df, times)

    result <- list(time=times, surv=t(survprobs), fullsurv=full_survprobs)
    attr(result, 'class') <- 'survfit.prev'
    result
}

print.survfit.prev <- function(object, ...) {
    cat("Survival probability calculated at", length(object$time), "timepoints, across", dim(object$surv)[1], "bootstraps.")
}

#' Summarises survival information at pre-specified years of interest on a
#' \code{survfit.prev} object.
#'
#' Survival probability is estimated as the mean of the bootstrapped survival curves at a specific
#' timepoint, with the 2.5% and 97.5% quantiles providing 95% confidence intervals. Survival probability
#' can only be estimated at timepoints less than the maximum survival time in the original fitting of
#' the \code{prevalence} object.
#'
#' @param object A \code{survfit.prev} object.
#' @param years A vector of years at which to estimate survival probability from the bootstrapped
#' survival curves.
#' @return None, instead displays the survival probabilities to screen as a side-effect.
#' @examples
#' data(prevsim)
#'
#' prev_obj <- prevalence(Surv(time, status) ~ age(age) + sex(sex) + entry(entrydate) + event(eventdate),
#'                        data=prevsim, num_years_to_estimate = c(5, 10), population_size=1e6,
#'                        start = "2005-09-01",
#'                        num_reg_years = 8, cure = 5)
#'
#' survobj <- survfit(prev_obj, newdata=list(age=65, sex=0))
#'
#' summary(survobj)
#'
#' summary(survobj, years=c(1, 3, 5, 7))
#'
#' @export summary.survfit.prev
summary.survfit.prev <- function(object, years=c(1, 3, 5), ...) {

    # Truncate years to the maximum number allowed
    max_years <- floor(dim(object$surv)[2] / 365.25)
    if (sum(years > max_years) > 0)
        message("Cannot estimate survival probabilities beyond the ", max_years, " years provided in the dataset.")
    years <- years[years <= max_years]
    days <- sapply(years, function(x) floor(x * 365.25))

    probs <- colMeans(as.matrix(object$surv[, days]))
    lower <- apply(object$surv[, days], 2, function(x) quantile(x, 0.025))
    upper <- apply(object$surv[, days], 2, function(x) quantile(x, 0.975))

    cat("Survival probability estimated using", dim(object$surv)[1], "bootstrap survival curves:\n")
    f <- mapply(function(x, y, l, u) cat(x, " year survival: ", y, " (", l, " - ", u, ") \n", sep=''),
                years, round(probs, 3), round(lower, 3), round(upper, 3))
}

.calculate_survival_probability <- function(coef, data, times) {
    scale <- exp(coef[-length(coef)] %*% data) + 0.000001  # Hack to ensure scale != 0
    shape <- 1 / coef[length(coef)]
    1 - pweibull(times, scale=scale, shape=shape)
}