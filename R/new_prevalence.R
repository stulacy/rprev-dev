library(flexsurv)
library(lubridate)
library(dplyr)
library(data.table)

INCIDENCE_MARGIN <- 1.5

# TODO Current hardcoded values:
#   - Column names
#   - Weibull distribution for the survival model (not going to touch custom incident and survival objects although in theory shouldn't
#       be too much work)
#   - Covariates used in survival modelling! Currently hardcoded as age and sex (could have these passed in as a variable? surv_covars,
#       or a formula? have separated formulae for incidence and survival?)
#   - Also hardcode age to force user death at 100

new_prevalence <- function(index, num_years_to_estimate,
                           data,
                           registry_start_date=NULL,
                           N_boot=1000,
                           population_size=NULL, proportion=100e3,
                           level=0.95,
                           precision=2, n_cores=1) {

    # This argument allows the user to specify when their registry started. I.e. it could have
    # started a month before received first incident case, in which case would want that taking into account when
    # estimating incidence rate and prevalence
    if (is.null(registry_start_date)) {
        registry_start_date <- min(data$entrydate)
    }

    index <- lubridate::ymd(index)
    registry_start_date <- lubridate::ymd(registry_start_date)
    sim_start_date <- index - lubridate::years(max(num_years_to_estimate))

    # If need simulation
    if (sim_start_date < registry_start_date) {
        sim_time <- difftime(index, sim_start_date, units='days')
        prev_sim <- new_sim_prevalance(data, index, sim_time, starting_date=sim_start_date, nsims=N_boot)
        # Create column indicating whether contributed to prevalence for each year of interest
        # For each year in Nyear TODO **THAT IS GREATER THAN NUMBER OF YEARS AVAILABLE IN REGISTRY**:
        for (year in num_years_to_estimate) {
            # Determine starting incident date
            starting_incident_date <- index - lubridate::years(year)

            # We'll create a new column to hold a binary indicator of whether that observation contributes to prevalence
            col_name <- paste0("prev_", year, "yr")

            # Determine prevalence as incident date is in range and alive at index date
            prev_sim[, (col_name) := as.numeric(incident_date > starting_incident_date & death_date > index)]
        }


    } else {
        prev_sim <- NULL
    }

    # Determine point estimates of prevalence by combining simulated and counted values
    names <- sapply(num_years_to_estimate, function(x) paste('y', x, sep=''))
    estimates <- lapply(setNames(num_years_to_estimate, names),
                        new_point_estimate,  # Function
                        prev_sim, index, data,
                        registry_start_date,
                        population_size, proportion, level, precision)


    # Return object
    object <- list(estimates=estimates, simulated=prev_sim,
                   counted=new_counted_prevalence(index, data, registry_start_date),
                   index_date=index,
                   proportion=proportion,
                   registry_start=registry_start_date)

    # Calculate covariate means and save
    # TODO What's this needed for?
    #mean_df <- data[, c(age_var, sex_var)]
    #mean_df <- apply(mean_df, 2, as.numeric)
    #object$means <- colMeans(mean_df)
    #object$y <- survobj

    if (!is.null(prev_sim)) {
        object$pval <- new_test_prevalence_fit(object)
    }

    attr(object, 'class') <- 'prevalence'
    object

}

new_point_estimate <- function(year, sim_results, index, registry_data, registry_start_date, population_size=NULL, proportion=100e3,
                               level=0.95, precision=2) {

    # See if need simulation if have less registry data than required
    initial_date <- index - years(year)
    need_simulation <- initial_date < registry_start_date

    #need_simulation <- (as.numeric(difftime(index, registry_start_date, units='weeks')) / 52.25) < year
    count_prev <- new_counted_prevalence(index, registry_data, max(initial_date, registry_start_date))

    # Estimate absolute prevalence as sum of simulated and counted
    if (need_simulation) {
        # Determine column name for this year
        col_name <- paste0("prev_", year, "yr")
        sim_contributions <- sim_results[incident_date < registry_start_date][, sum(get(col_name)), by=sim][[2]]  # Return results column
        the_estimate <- count_prev + mean(sim_contributions)
    } else {
        the_estimate <- count_prev
    }

    result <- list(absolute.prevalence=the_estimate)

    if (!is.null(population_size)) {
        raw_proportion <- the_estimate / population_size
        the_proportion <- proportion * raw_proportion

        # TODO Have some way of testing whether had simulated data or not
        if (! need_simulation) {
            se <- (raw_proportion * (1 - raw_proportion)) / population_size
        } else {

            # Counted prevalence SE
            raw_proportion_n <- count_prev / population_size
            std_err_1 <- sqrt((raw_proportion_n * (1 - raw_proportion_n)) / population_size)

            # Simulated prevalence SE
            sim_prev_sd <- sd(sim_contributions)
            std_err_2 <- sim_prev_sd / population_size

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



# Requires columns with entrydate, eventdate, status
# Start date allows users to specify how long estimating prevalence for, as otherwise including
# all contributions in data set
# TODO Make start date optional like in prevalence
new_counted_prevalence <- function(index, data, start_date) {
    data %>%
        filter(entrydate >= start_date, entrydate < index) %>%
        mutate(dead_at_index = ifelse(eventdate > index, 0, status)) %>%
        summarise(sum(1 - dead_at_index)) %>%     # Prevalence is number of those still alive or censored at index
        unname() %>%
        as.numeric()
}

new_sim_prevalance <- function(data, index, number_incident_days, starting_date=NULL, inc_model=NULL, surv_model=NULL, nsims=1000,
                               dist='weibull')  # Currently does nothing
{
    full_data <- data

    # If inc_model is null then make one ourselves, must be a func of time and N and returns interarrival times
    if (is.null(inc_model)) {
        # TODO Make the incidence model accept the registry starting date as an argument too
        inc_model <- fit_exponential_incidence(data)
    }

    # If surv_model is null then make one ourselves, must be a function of data frame with age and sex values that returns draws from
    # survival distribution
    # TODO Allow for the survival distribution to be chosen as one that ISN'T 3-parameter
    if (is.null(surv_model)) {
        surv_model <- flexsurvreg(Surv(time, status) ~ age + sex, data=data, dist='weibull')
    }

    all_results <- replicate(nsims, {
        # bootstrap dataset
        data <- full_data[sample(seq(nrow(full_data)), replace=T), ]

        # fit incidence and survival models.
        bs_inc <- eval(inc_model$call)
        bs_surv <- eval(surv_model$call)

        # Estimate how many people will be incident in this time frame
        initial_num_inds <- INCIDENCE_MARGIN * expected_incidence(bs_inc, number_incident_days)

        # draw incidence times (time in days after the starting date)
        entry_times <- draw_interarrival_time(bs_inc, initial_num_inds)
        entry_times <- entry_times[entry_times < number_incident_days]
        num_inds <- length(entry_times)

        # draw covariate values from model (TODO maybe require that model has a "terms" attribute?)
        # TODO Expand this to determine the covariates in the model rather than being hardcoded
        newdata <- data.frame(age=sample(full_data$age, num_inds, replace=T),
                              sex=sample(full_data$sex, num_inds, replace=T)
                              )

        # draw event times
        event_times <- draw_time_to_death.flexsurvreg(surv_model, newdata)

        # Turn into Data Table (far quicker than data frame)
        dt <- data.table(newdata)
        dt[, c("time_to_entry", "time_to_death") := list(entry_times, event_times)]
        dt
    }, simplify=FALSE)

    # Combine into single table
    results <- rbindlist(all_results, idcol='sim')

    # truncate death at age 100
    results[(age*365.25 + time_to_death) > 36525, time_to_death := 36525 - age*365.25]

    # Add the data into the equation, i.e. incidence date, death date
    if (!is.null(starting_date)) {
        results[, c('incident_date', 'death_date') := list(as.Date(starting_date + time_to_entry),
                                                           as.Date(starting_date + time_to_entry + time_to_death)
                                                           )]
    }

    # Drop intermediary columns
    results
}

# TODO Put these in a survival script somewhere
draw_time_to_death <- function(object, newdata, ...) UseMethod("draw_time_to_death")

draw_time_to_death.flexsurvreg <- function(object, newdata) {
    # Obtain location parameter (i.e. one that covariates act on)
    loc_param <- object$dlist$location

    # Obtain linear predictor for the location parameter
    # Obtain name of covariates
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
