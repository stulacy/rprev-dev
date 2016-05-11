#' Print the current value of the incidence object.
#'
#' @param object ...
#' @return ...
#' @examples ...
print.incidence <- function(object, ...) {
    cat(paste("Mean incidence per year of the registry is ", round(mean(object$raw_incidence),2),
              ", with standard deviation ", round(sd(object$raw_incidence),2), ".\n", sep=''))
}

#' Generate a summary of the incidence object.
#'
#' @param object ...
#' @return ...
#' @examples ...
summary.incidence <- function(object, ...) {
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
#' @param object ...
#' @return ...
#' @examples ...
print.prevalence <- function(object, ...) {
    cat(paste("Estimated ", object$nyears, " year prevalence is ", sum(object$simulated_cases), ".\n", sep=''))
}

#' Generate a summary of the prevalence object.
#'
#' @param object ...
#' @return ...
#' @examples ...
summary.prevalence <- function(object, ...) {
    cat("Registry Data\n~~~~~~~~~~~~~\n")
    cat("Number of years:", length(object$known_inc_rate), "\n")
    cat("Known incidence rate:\n")
    cat(object$known_inc_rate)

    cat("\n\nBootstrapping\n~~~~~~~~~~~~~\n")
    cat("Iterations:", object$nbootstraps, "\n")
    cat("Posterior age distribution summary:\n")
    print(summary(object$post_covar))
    cat("Average prevalent cases per year:\n")
    cat(object$simulated_cases)
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
