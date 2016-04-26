print.prevalence <- function(object, ...) {
    cat(paste("Estimated ", object$nyears, " year prevalence is ", sum(object$cases_avg), ".\n", sep=''))
}

summary.prevalence <- function(object, ...) {
    cat("Registry Data\n~~~~~~~~~~~~~\n")
    cat("Number of years:", length(object$known_inc_rate), "\n")
    cat("Known incidence rate:\n")
    cat(object$known_inc_rate)
    
    cat("\n\nBootstrapping\n~~~~~~~~~~~~~\n")
    cat("Iterations:", object$nbootstraps, "\n")
    cat("Posterior age distribution summary:\n")
    print(summary(object$post))
    cat("Average prevalent cases per year:\n")
    cat(object$cases_avg)
}