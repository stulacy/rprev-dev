#' Summarise disease incidence.
#'
#' Calculates incidence by year of the registry data, along with mean incidence
#' with confidence intervals. A smoothed cumulative incidence function is fit to
#' the data for inspecting deviations in the registry data from a homogeneous
#' Poisson process.
#'
#' @inheritParams prevalence
#' @param entry Vector of diagnosis dates for each patient in the registry in
#'   the format YYYY-MM-DD.
#' @param start Date from which incident cases are included in the format
#'   YYYY-MM-DD. Defaults to the earliest entry date.
#' @param num_reg_years The number of years of the registry for which incidence
#'   is to be calculated. Defaults to using all available complete years.
#' @param df The desired degrees of freedom of the smooth.
#' @param precision The number of decimal places required.
#' @param level The desired confidence interval width.
#' @return An S3 object of class \code{incidence} with the following attributes:
#'   \item{raw_incidence}{Vector of absolute incidence values for each included
#'   year of the registry, as generated using \code{\link{yearly_incidence}}.}
#'   \item{ordered_diagnoses}{Vector of times (days) between diagnosis date and
#'   the earliest date of inclusion in the registry, ordered shortest to
#'   longest.} \item{smooth}{Smooth fitted to the cumulative incidence data.}
#'   \item{index_dates}{Dates delimiting the years in which incidence is
#'   calculated.} \item{mean}{List containing mean incidence per 100K with
#'   confidence intervals. See \link{mean_incidence_rate}.} \item{dof}{Degrees
#'   of freedom of the smooth.}
#' @examples
#' data(prevsim)
#'
#' \dontrun{
#' incidence(prevsim$entrydate, 1e6)
#'
#' incidence(prevsim$entrydate, 1e6, start = "2004-01-30", num_reg_years = 9)
#' }
#'
#' @export
#' @family incidence functions
incidence <- function(entry, population_size, start=NULL, num_reg_years=NULL,
                      df=6, precision=2, level=0.95){

    if (is.null(start))
        start <- min(entry)

    if (is.null(num_reg_years))
        num_reg_years <- floor(as.numeric(difftime(max(entry), start) / DAYS_IN_YEAR))

    reg_years <- determine_yearly_endpoints(start, num_reg_years)
    index_date <- reg_years[length(reg_years)]

    entry_trunc <- entry[entry >= start & entry < index_date]

    # Slightly confused that the following are not all integers:
    diags <- sort(as.numeric(difftime(entry_trunc, min(entry_trunc), units='days')))
    smo <- smooth.spline(diags, seq(length(diags)), df=df)
    raw_inc <- yearly_incidence(entry, start_date=start, num_years=num_reg_years)

    object <- list(raw_incidence=raw_inc,
                   ordered_diagnoses=diags,
                   smooth = smo,
                   index_dates = reg_years,
                   mean=mean_incidence_rate(raw_inc, population_size=population_size,
                                            precision=precision, level=level),
                   pvals=test_incidence_fit(raw_inc),
                   dof=df)
    attr(object, 'class') <- 'incidence'
    object
}

#' DEPRECATED: Please use \code{\link{yearly_incidence}} instead.
#'
#' Disease incidence.
#'
#' Calculates yearly incidence for the available registry data.
#'
#' @inheritParams incidence
#' @return Vector of length num_reg_years of integers, representing the number
#'   of absolute incidence values for each included year of the registry.
#' @examples
#' data(prevsim)
#'
#' @export
#' @family incidence functions
#' @seealso \link{yearly_incidence}
raw_incidence <- function(entry, start=NULL, num_reg_years=NULL) {
    .Deprecated("yearly_incidence")
    yearly_incidence(entry, start_date=start, num_years=num_reg_years)

}

#' Disease incidence.
#'
#' Calculates yearly incidence for the available registry data.
#'
#' @inheritParams incidence
#' @param start_date The initial date in the \code{entry} vector to start estimating incidence from.
#' @param num_years The number of complete years to calculate incidence over. Defaults to the number of complete
#' years of registry data available in \code{entry}.
#' @param end_date The ending date in the \code{entry} vector to estimate incidence counting back from.
#' If both \code{end_date} and \code{start_date} are specified then \code{start_date} takes precedence.
#' @return Vector of length \code{num_years} of integers, representing the number
#'   of absolute incidence values for each included year of the registry.
#' @examples
#' data(prevsim)
#'
#' yearly_incidence(prevsim$entrydate, start_date="2004-01-01", num_years=8)
#' yearly_incidence(prevsim$entrydate)
#' yearly_incidence(prevsim$entrydate, start_date="2005-05-01", num_years=5)
#' yearly_incidence(prevsim$entrydate, start_date="2005-05-01")
#' yearly_incidence(prevsim$entrydate, num_years=5, end_date="2015-05-01")
#'
#' @export
#' @family incidence functions
yearly_incidence <- function(entry, start_date=NULL, num_years=NULL, end_date=NULL) {

    if (!is.null(start_date)) { # Having the start date takes priority
        if (is.null(num_years)) {
            if (is.null(end_date)) {
                end_date <- max(entry)
            }
            num_years <- floor(as.numeric(difftime(max(entry), start_date) / DAYS_IN_YEAR))
        }
        registry_years <- determine_yearly_endpoints(start_date, num_years)
    } else {
        if (!is.null(end_date)) { # Having end date takes second priority
            if (is.null(num_years)) {
                num_years <- floor(as.numeric(difftime(end_date, min(entry)) / DAYS_IN_YEAR))
            }
            registry_years <- determine_yearly_endpoints(end_date, num_years, direction='backwards')
        } else { # No start date, no end date
            start_date <- min(entry)
            if (is.null(num_years)) {
                num_years <- floor(as.numeric(difftime(max(entry), start_date) / DAYS_IN_YEAR))
            }
            registry_years <- determine_yearly_endpoints(start_date, num_years)
        }
    }

    # Force date in case data isn't supplied correctly
    entry <- as.Date(entry)
    per_year <- vapply(seq(num_years),
                       function(i) sum(entry >= registry_years[i] & entry < registry_years[i+1]),
                       integer(1))

    if(sum(per_year) < 30) warning("Warning: low number of incident cases.")
    per_year
}

#' Mean disease incidence.
#'
#' Calculates the average incidence rate per one hundred thousand with
#' confidence intervals for the given registry data.
#'
#' @inheritParams incidence
#' @param raw_inc Vector of incidence values by year, as returned by
#'   \code{\link{yearly_incidence}}.
#' @return A list with the following values:
#'
#'   \item{absolute}{Overall incidence for the period of interest as a single
#'   double} \item{per100K}{Incidence for the period of interest per one hundred
#'   thousand} \item{per100K.lower}{Lower bounds of the specified confidence
#'   level on the per one hundred thousand estimate} \item{per100K.upper}{Upper
#'   bounds of the specified confidence level on the per one hundred thousand
#'   estimate}
#'
#' @examples
#' data(prevsim)
#'
#' rawinc <- yearly_incidence(prevsim$entrydate)
#' mean_incidence_rate(rawinc, population_size=3.5e6)
#'
#' rawinc2 <- yearly_incidence(prevsim$entrydate, start_date="2005-05-01", num_years=5)
#' mean_incidence_rate(rawinc2, population_size=3.5e6)
#'
#' @export
#' @family incidence functions
mean_incidence_rate <- function(raw_inc, population_size, precision = 2, level=0.95){

    mean_rate <- mean(raw_inc)
    z_conf <- qnorm((1+level)/2)

    CI <- 100e3 * (z_conf * sqrt(mean_rate) / length(raw_inc)) / population_size
    est <- 100e3 * mean_rate / population_size

    object <- list(absolute=mean_rate, per100K=est, per100K.lower=est-CI, per100K.upper=est+CI)
    lapply(object, round, precision)
}


#' Visualise disease incidence.
#'
#' Plots a comparison between the smoothed daily incidence function and actual
#' incidence.
#'
#' This function generates a plot from the cumulative incidence object. The
#' incidence rate per year of the registry is shown in red. Mean incidence rate
#' is shown as a solid blue line, with the confidence interval shown in dashed
#' blue lines. The smooth fitted to the cumulative incidence data is shown in
#' green.
#' @param x An \code{incidence} object.
#' @param level The desired confidence interval width.
#' @param ... Arguments passed to \code{plot}.
#' @return An object of class \code{ggplot}.
#' @examples
#' data(prevsim)
#'
#' \dontrun{
#' inc <- incidence(prevsim$entrydate, population_size=1e6, start = "2004-01-30", num_reg_years = 9)
#'
#' plot(inc)
#' }
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 scale_colour_manual
#' @export
#' @family incidence functions
plot.incidence <- function(x, level=0.95, ...){
    raw_incidence <- x$raw_incidence
    mean_rate <- mean(raw_incidence)
    day_mean_rate <- mean_rate / DAYS_IN_YEAR

    z_conf <- qnorm((1+level)/2)
    CI_lim <- z_conf * sqrt(mean_rate)/DAYS_IN_YEAR
    num_reg_years <- length(raw_incidence)

    inc_rate <- data.frame(inc=raw_incidence/DAYS_IN_YEAR, day=as.Date(x$index_dates[-length(x$index_dates)]) + 182, col='r')
    pred_rate <- predict(x$smooth, seq(num_reg_years*DAYS_IN_YEAR), deriv=1)
    smooth_rate <- data.frame(rate=pred_rate$y, day=as.Date(x$index_dates[1]) + pred_rate$x, col='g')
    mean_rate <- data.frame(mean=day_mean_rate, upper=day_mean_rate+CI_lim,
                            lower=day_mean_rate-CI_lim, col='b')
    ci_diff <- 0.5 * (mean_rate$upper - mean_rate$lower)

    p <- ggplot2::ggplot() +
            ggplot2::geom_point(data=inc_rate, ggplot2::aes_string(x='day', y='inc', colour='col')) +
            ggplot2::geom_line(data=inc_rate, ggplot2::aes_string(x='day', y='inc', colour='col'), size=1) +
            ggplot2::geom_line(data=smooth_rate, ggplot2::aes_string(x='day', y='rate', colour='col'),  size=1) +
            ggplot2::geom_hline(data=mean_rate, ggplot2::aes_string(yintercept='mean', colour='col')) +
            ggplot2::geom_hline(data=mean_rate, ggplot2::aes_string(yintercept='upper', colour='col'),
                                linetype='dashed') +
            ggplot2::geom_hline(data=mean_rate, ggplot2::aes_string(yintercept='lower', colour='col'),
                                linetype='dashed') +
            ggplot2::labs(x='Year', y='Daily incidence rate') +
            ggplot2::theme_bw() +
            ggplot2::theme(legend.position='bottom') +
            ggplot2::scale_colour_manual(name='Data',
                                         values=c('r'='red', 'g'='green', 'b'='#0080ff'),
                                         breaks=c('r', 'g', 'b'),
                                         labels=c('Actual incidence', 'Smoothed incidence', 'Mean actual incidence'))
    p
}

#' Visualise incidence Poisson assumption.
#'
#' Plots a figure detailing the deviation of incidence rate from an homogeneous
#' Poisson process using simulation.
#'
#' This function generates a plot where the red line is the deviation of
#' cumulative diagnoses from the fitted smooth taken from the incidence object.
#' In order to ascertain if the deviation is within reasonable bounds,
#' simulation is used. The grey lines represent N draws from a uniform
#' distribution between 0 and the last day of diagnosis, where N is the number
#' of incident cases, plotted as deviations from a smooth fitted to them. The
#' blue lines are 95\% confidence intervals drawn from the simulated data.
#'
#' @param object An \code{incidence} object.
#' @param N_sim Number of draws from a homogeneous Poisson process.
#' @param level The desired confidence interval width.
#' @param samples_per_bin Number of samples per bin.
#' @param max_bins Maximum number of bins.
#' @return An object of class \code{ggplot}.
#' @examples
#' data(prevsim)
#'
#' \dontrun{
#' inc <- incidence(prevsim$entrydate, 1e6, start = "2004-01-30",
#'                  num_reg_years = 9)
#'
#' plot_incidence_fit(inc)
#' }
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 geom_ribbon
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 aes_string
.plot_incidence_fit <- function(object, N_sim=1000, level=0.95, samples_per_bin=10, max_bins=200){
    diags <- object$ordered_diagnoses
    N <- length(diags)

    boot_out <- lapply(seq(N_sim), function(i) {
        x <- sort(runif(N, 0, max(diags)))
        smo <- smooth.spline(x, 1:N, df=object$dof)
        data.frame(y=seq(N) - predict(smo, x)$y, x=x)
    })

    bootstraps = do.call(rbind, boot_out)

    num_bins_init <- floor(N / samples_per_bin)
    num_bins <- if (num_bins_init > max_bins) max_bins else num_bins_init
    bin_segments <- seq(0, max(diags), by=max(diags) / num_bins)
    bin_segments[length(bin_segments)] <- bin_segments[length(bin_segments)] + 1  # Expand last bin to allow for max data point

    bin_cis = lapply(seq(num_bins-1), function(i) {
        vals <- lapply(boot_out, function(bootstrap) {
            bootstrap[bootstrap$x >= bin_segments[i] & bootstrap$x < bin_segments[i+1], 'y']
        })
        avgs <- sapply(vals, mean)
        data.frame(x=(sum(bin_segments[c(i, i+1)])) / 2, upper=quantile(avgs, (1+level)/2, na.rm=T),
                   lower=quantile(avgs, 1-(1+level)/2, na.rm=T))
    })
    bin_cis_df <- do.call(rbind, bin_cis)

    p <- ggplot2::ggplot() +
        ggplot2::geom_vline(data=bin_cis_df, ggplot2::aes_string(xintercept='x'), alpha=0.2, linetype='dotted') +
        ggplot2::geom_ribbon(data=bin_cis_df, ggplot2::aes_string(x='x', ymin='lower', ymax='upper'), alpha=0.1) +
        ggplot2::theme_bw() +
        ggplot2::geom_line(data=data.frame(x=diags, y=seq(N)-predict(object$smooth, diags)$y),
                           ggplot2::aes_string(x='x', y='y'),
                           colour='orange', size=1) +
        ggplot2::labs(x="Days", y="Deviation from smooth")
    p
}

#' @export
print.incidence <- function(x, ...) {
    cat("Cumulative incidence object with", length(x$raw_incidence), "years of data.\n")
    cat("Smooth fitted using", x$dof, "degrees of freedom.\n")

}

#' @export
summary.incidence <- function(object, ...) {
    cat("Number of years of registry data:", length(object$raw_incidence), "\n")

    cat("\nIncidence\n~~~~~~~~~\n")
    cat("Known incidence by year:", object$raw_incidence, "\n")
    cat("Diagnoses (time since registry began):\n")
    print(summary(object$ordered_diagnoses))
    cat("p-values for over/under dispersion:", object$pvals, "\n")

    cat("\nFitted smooth:\n~~~~~~~~~~~~~\n")
    print(object$smooth)
}
