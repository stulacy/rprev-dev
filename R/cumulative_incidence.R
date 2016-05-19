#' Generates a smoothed cumulative incidence function, useful for inspecting deviations in the registry data
#' and compare with the cumulative incidence for constant diagnosis rate.
#'
#' @param entry Vector of diagnosis dates for each patient in the registry in the format YYYY-MM-DD.
#' @param start Date from which incident cases are included in the format YYYY-MM-DD,
#' defaults to the earliest entry date.
#' @param num_reg_years The number of years of the registry for which incidence is to be calculated,
#' defaults to using all available complete years.
#' @param df The desired degrees of freedom for the smoothening function.
#' @return An S3 object of class \code{cincidence} with the following attributes:
#' \item{raw_incidence}{Vector of absolute incidence values for each included year of the registry,
#' as generated using \code{incidence()}.}
#' \item{ordered_diagnoses}{Vector of times (days) between diagnosis date and the earliest date of
#' inclusion in the registry, ordered shortest to longest.}
#' \item{smooth}{Smooth fitted to the cumulative incidence data.}
#' \item{index_dates}{Dates delimiting the years in which incidence is calculated.}
#' @examples
#' data(prevsim)
#'
#' c_inc <- cumulative_incidence(previm$entrydate, start = "2004-01-30", num_reg_years = 9)
#' ordered_diagnoses <- c_inc$ordered_diagnoses
#'
#' # Plot the smooth against the ordered diagnoses to see the model fit
#' plot(ordered_diagnoses, c_inc$cumulative_incidence, pch=20,
#' cex=0.7, xlab="days", ylab="cumulative diagnoses")
#' abline(a=0, b=length(ordered_diagnoses)/ordered_diagnoses[length(ordered_diagnoses)],
#' col="red", lwd=2)
#' lines(c_inc$smooth, col="green", lwd=2)
#'
#' # The following plot shows the deviation of the raw data from the fitted smooth by day:
#' plot(ordered_diagnoses, seq(length(ordered_diagnoses)) - predict(c_inc$smooth, ordered_diagnoses)$y,
#' type="l", xlab="days", ylab="deviation from smooth")
#' @export cumulative_incidence
cumulative_incidence <- function(entry, start = NULL, num_reg_years = NULL, df=6){

    if (is.null(start))
        start <- min(entry)

    if (is.null(num_reg_years))
        num_reg_years <- floor(as.numeric(difftime(max(entry), start) / 365.25))

    reg_years <- determine_registry_years(start, num_reg_years)
    index_date <- reg_years[length(reg_years)]

    entry_trunc <- entry[entry >= start & entry < index_date]

    # Slightly confused that the following are not all integers:
    diags <- sort(as.numeric(difftime(entry_trunc, min(entry_trunc), units='days')))
    smo <- smooth.spline(diags, seq(length(diags)), df=df)

    cumulative_inc_out <- list(raw_incidence = incidence(entry, start, num_reg_years),
                               ordered_diagnoses = diags,
                               smooth = smo,
                               index_dates = reg_years)
    attr(cumulative_inc_out, 'class') <- 'cincidence'
    cumulative_inc_out
}

#' Plots a comparison between the smoothed daily incidence function and the actual.
#'
#' This function generates a plot from the cumulative incidence object. The incidence rate per
#' year of the registry is shown in red. Mean incidence rate is shown as a solid blue line,
#' with the 95\% confidence interval shown with dashed blue lines. The smooth fitted to the
#' cumulative incidence data is shown in green.
#' @param object A \code{cincidence} object.
#' @param level The desired confidence interval width.
#' @return Plot of incidence rate, confidence interval and smoothed incidence function as a side-effect.
#' @examples
#' data(prevsim)
#'
#' c_inc <- cumulative_incidence(previm$entrydate, start = "2004-01-30", num_reg_years = 9)
#'
#' inspect_incidence(c_inc)
#' @export inspect_incidence
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 ylim
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 scale_colour_manual
inspect_incidence <- function(object, level=0.95){
    raw_incidence <- object$raw_incidence
    mean_rate <- mean(raw_incidence)
    day_mean_rate <- mean_rate / 365

    z_conf <- qnorm((1+level)/2)
    CI_lim <- z_conf * sqrt(mean_rate)/365
    num_reg_years <- length(raw_incidence)

    inc_rate <- data.frame(inc=raw_incidence/365, day=as.Date(object$index_dates[-length(object$index_dates)]) + 182)
    pred_rate <- predict(object$smooth, seq(num_reg_years*365), deriv=1)
    smooth_rate <- data.frame(rate=pred_rate$y, day=as.Date(object$index_dates[1]) + pred_rate$x)
    mean_rate <- data.frame(mean=day_mean_rate, upper=day_mean_rate+CI_lim,
                            lower=day_mean_rate-CI_lim)
    ci_diff <- 0.5 * (mean_rate$upper - mean_rate$lower)

    ggplot2::ggplot() +
        ggplot2::geom_point(data=inc_rate, ggplot2::aes(x=day, y=inc, colour='r')) +
        ggplot2::geom_line(data=inc_rate, ggplot2::aes(x=day, y=inc, colour='r'), size=1) +
        ggplot2::geom_line(data=smooth_rate, ggplot2::aes(x=day, y=rate, colour='g'),  size=1) +
        ggplot2::geom_hline(data=mean_rate, ggplot2::aes(yintercept=mean, colour='b')) +
        ggplot2::geom_hline(data=mean_rate, ggplot2::aes(yintercept=upper, colour='b'), linetype='dashed') +
        ggplot2::geom_hline(data=mean_rate, ggplot2::aes(yintercept=lower, colour='b'), linetype='dashed') +
        ggplot2::labs(x='Year', y='Daily incidence rate') +
        ggplot2::theme_bw() +
        ggplot2::ylim(mean_rate$lower-ci_diff, mean_rate$upper+ci_diff) +
        ggplot2::scale_colour_manual(name='Data',
                                     values=c('r'='red', 'g'='green', 'b'='#0080ff'),
                                     breaks=c('r', 'g', 'b'),
                                     labels=c('Actual incidence', 'Smoothed incidence', 'Mean actual incidence'))
}

#' Inspect consistency of incidence rate with an homogeneous Poisson process using simulation.
#'
#' This function generates a plot where the red line is the deviation of cumulative diagnoses
#' from the fitted smooth taken from the incidence object. In order to ascertain if the deviation
#' is within reasonable bounds, simulation is used. The grey lines represent N draws from a
#' uniform distribution between 0 and the last day of diagnosis, where N is the number of incident
#' cases, plotted as deviations from a smooth fitted to them. The blue lines are 95\% confidence
#' intervals drawn from the simulated data.
#'
#' @param object A \code{cincidence} object.
#' @param N_sim Number of draws from a homogeneous Poisson process.
#' @param level The desired confidence interval width.
#' @param df The desired degrees of freedom for the smoothening function.
#' @return Plot of the smoothed incidence function and corresponding deviations.
#' @examples
#' poisson_incidence_sim(c_inc)
poisson_incidence_sim <- function(object, N_sim=1000, level=0.95, df=4){

  diags <- object$ordered_diagnoses
  N <- length(diags)
  boot_out <- matrix(NA, nrow = N_sim, ncol = N)

  for (i in 1:N_sim){
    x <- sort(runif(N, 0, max(diags)))
    the_smo <- smooth.spline(x, 1:N, df=df)
    boot_out[i, ] <- (1:N) - predict(the_smo, x)$y
  }

  plot(NA, xlim=c(0, max(diags)), ylim=c(-0.8*max(boot_out),0.8*max(boot_out)), xlab="days", ylab="deviation from smooth")
  sapply(seq(N_sim),
         function(i) lines(x, boot_out[i,], col="grey"))

  lines(diags, seq(length(diags)) - predict(object$smooth, diags)$y, col="red")
  lines(x, apply(boot_out, 2, quantile, probs=(1+level)/2), col="blue")
  lines(x, apply(boot_out, 2, quantile, probs=1-((1+level)/2)), col="blue")

}

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
