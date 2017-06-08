library(flexsurv)
library(data.table)

INCIDENCE_MARGIN <- 1.5

# TODO Current hardcoded values:
#   - Column names
#   - Weibull distribution for the survival model (not going to touch custom incident and survival objects although in theory shouldn't
#       be too much work)
#   - Covariates used in survival modelling! Currently hardcoded as age and sex (could have these passed in as a variable? surv_covars,
#       or a formula? have separated formulae for incidence and survival?)
#   - Also hardcode age to force user death at 100

new_sim_prevalance <- function(data, index='2015-08-31', Nyears=c(10), inc_model=NULL, surv_model=NULL, nsims=1000,
                               dist='weibull')  # Currently does nothing
{
    full_data <- data
    ndays <- max(Nyears * 365.25)
    starting_date <- as.Date(index) - Nyears * 365.25

    # TODO Somewhere sort out how to calculate multiple starting_dates, so can obtain prevalence for each N year required

    # If inc_model is null then make one ourselves, must be a func of time and N and returns interarrival times
    if (is.null(inc_model)) {
        inc_model <- fit_exponential_incidence(data)
    }

    # If surv_model is null then make one ourselves, must be a function of data frame with age and sex values that returns draws from
    # survival distribution
    # TODO Allow for the survival distribution to be chosen as one that ISN'T 3-parameter
    if (is.null(surv_model)) {
        surv_model <- flexsurvreg(Surv(time, status) ~ age + sex, data=data, dist='weibull')
    }

    all_results <- lapply(1:nsims, function(i) {
        # bootstrap dataset
        data <- full_data[sample(seq(nrow(full_data)), replace=T), ]

        # fit incidence and survival models
        bs_inc <- eval(inc_model$call)
        bs_surv <- eval(surv_model$call)

        # Estimate how many people will be incident in this time frame
        initial_num_inds <- INCIDENCE_MARGIN * expected_incidence(bs_inc, ndays)

        # draw incidence times (time in days after the starting date (index date - Nyears))
        entry_times <- draw_interarrival_time(bs_inc, initial_num_inds)
        entry_times <- entry_times[entry_times < ndays]
        num_inds <- length(entry_times)

        # draw age and sex values
        # TODO Expand this to determine the covariates in the model rather than being hardcoded
        newdata <- data.frame(age=sample(full_data$age, num_inds, replace=T),
                              sex=sample(full_data$sex, num_inds, replace=T)
                              )

        # draw event times
        event_times <- draw_time_to_death.flexsurvreg(surv_model, newdata)

        # Turn into Data Table (far quicker than data frame)
        dt <- data.table(newdata)
        dt[, c("time_to_entry", "time_to_death") := list(entry_times, event_times)]
        # truncate death at age 100
        dt[(age*365.25 + time_to_death) > 36525, time_to_death := 36525 - age*365.25]
        dt
    })

    # Combine into single table
    results <- rbindlist(all_results, idcol='sim')

    # Add the data into the equation, i.e. incidence date, death date
    results[, c('incident_date', 'death_date') := list(starting_date + time_to_entry,
                                                       starting_date + time_to_entry + time_to_death
                                                       )]

    # Create column for prevalent at index
    # TODO Have this iterated (or grouped by) each year of prevalence we're interested in
    results[, prevalent := as.numeric(death_date > index)]

    # Can view the number of prevalent cases by each simulation
    prev_summary <- results[,. (num_prev = sum(prevalent)), by=sim]

    # And also the mean and CIs
    mean(prev_summary$num_prev)
    quantile(prev_summary$num_prev, c(0.025, 0.975))
}

# TODO Put these in a survival script somewhere
draw_time_to_death <- function(object, newdata, ...) UseMethod("draw_time_to_death")

draw_time_to_death.flexsurvreg <- function(object, newdata) {
    # Obtain location parameter (i.e. one that covariates act on)
    loc_param <- object$dlist$location

    # Obtain linear predictor for the location parameter
    # Obtain name of covariates (typically just 'age' and 'sex')
    covars <- colnames(newdata)
    formula <- as.formula(paste("~", paste(covars, collapse='+')))
    # Expand categorical variables in newdata frame
    wide_df <- model.matrix(formula, newdata)
    # Cross-multiply by coefficients to obtain linear predictor for location parameter
    lps <- wide_df %*% object$coefficients[c(loc_param, colnames(wide_df)[2:ncol(wide_df)])]

    # Inverse transform params to obtain them the way they are specified in R
    loc_value <- object$dlist$inv.transforms

    # This obtains the index which the location parameter is referred to in the list of params
    loc_index <- which(object$dlist$pars == object$dlist$location)
    scale_index <- which(object$dlist$pars != object$dlist$location)

    # Calculate the values with the inverse transforms for each of these parameters and
    # therefore also the time-to-death
    loc_vals <- object$dlist$inv.transforms[[loc_index]](lps)

    if (!is.null(scale_index)) {
        scale_val <- object$dlist$inv.transforms[[scale_index]](object$coefficients[scale_index])
        if (scale_index < loc_index) {
            object$dfns$r(nrow(newdata), scale_val, loc_vals)
        } else {
            object$dfns$r(nrow(newdata), loc_vals, scale_val)
        }
    } else {
        object$dfns$r(nrow(newdata), loc_vals)
    }

}

# TODO Put these in an incidence script somewhere
fit_exponential_incidence <- function(data) {
    rate <- length(data$entrydate) / as.numeric(max(as.Date(data$entrydate)) - min(as.Date(data$entrydate)))
    obj <- list(rate=rate, call=match.call())
    class(obj) <- "expinc"
    obj
}

draw_interarrival_time <- function(object, N, ...) UseMethod("draw_interarrival_time")

draw_interarrival_time.expinc <- function(object, N) {
    cumsum(rexp(N, object$rate))
}

expected_incidence <- function(object, timeframe, ...) UseMethod("expected_incidence")

expected_incidence.expinc <- function(object, timeframe) {
    timeframe * object$rate
}
