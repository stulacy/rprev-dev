#' Generate a summary of the cumulative_incidence object.
#'
#' @param object A \code{cumulative_incidence} object.
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

#' Print the current value of the prevalence object.
#'
#' @param object A \code{prevalence} object.
print.prevalence <- function(object, ...) {
    cat("Estimated prevalence per", object$proportion, "at", object$index_date, "\n")
    lapply(names(object$estimates), function(x) {
        year <- strsplit(x, 'y')[[1]][2]
        prev_est <- object$estimates[[x]][2]
        cat(paste(year, "years:", prev_est, "\n"))
    })

   # cat(paste("Estimated ", object$nyears, " year prevalence is ", sum(object$mean_yearly_contributions), ".\n", sep=''))
}

#' Generate a summary of the prevalence object.
#'
#' @param object A \code{prevalence} object.
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
    cat(round(rev(object$simulated$mean_yearly_contributions)))
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
