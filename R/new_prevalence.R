#library(flexsurv)
#library(lubridate)
#library(dplyr)
#library(data.table)

INCIDENCE_MARGIN <- 1.5
MIN_INCIDENCE <- 10

#' @import data.table
#' @export
new_prevalence <- function(index, num_years_to_estimate,
                           data,
                           inc_formula=NULL,
                           inc_model=NULL,
                           surv_formula=NULL,
                           surv_model=NULL,
                           registry_start_date=NULL,
                           death_column=NULL,
                           incident_column=NULL,
                           age_column='age',
                           status_column='status',
                           N_boot=1000,
                           population_size=NULL, proportion=100e3,
                           level=0.95,
                           dist='weibull',
                           precision=2, n_cores=1) {

    if (is.null(incident_column)) {
        if (!is.null(inc_formula)) {
            incident_column <- all.vars(update(inc_formula, .~0))
        }
    }

    # Is it right for surv_formula to have precedence over surv_formula?
    if (!is.null(surv_formula)) {
        surv_LHS <- all.vars(update(surv_formula, .~0))
        status_column <- surv_LHS[length(surv_LHS)]
    }

    # Form formula for counted prevalence
    # extract entry column from incidence formula
    if (is.null(death_column)) {
        message("death_column not provided so prevalence cannot be counted over the registry. Estimates will be solely from simulation.")
        counted_formula <- NULL
    } else if (is.null(incident_column)) {
        message("incident_column not provided so prevalence cannot be counted over the registry. Estimates will be solely from simulation.")
        counted_formula <- NULL
    } else {
        counted_formula <- as.formula(paste(death_column, incident_column, sep='~'))
    }

    if (!is.null(incident_column)) {
        data[[incident_column]] <- lubridate::ymd(data[[incident_column]])
    }
    if (!is.null(death_column)) {
        data[[death_column]] <- lubridate::ymd(data[[death_column]])
    }

    # This argument allows the user to specify when their registry started. I.e. it could have
    # started a month before received first incident case, in which case would want that taking into account when
    # estimating incidence rate and prevalence
    if (is.null(registry_start_date)) {
        if (is.null(incident_column)) {
            stop("Error: Unknown registry starting date. Please either provide 'registry_start_date' or the incident column in 'inc_formula'.")
        }
        registry_start_date <- min(data[[incident_column]])
    }

    index <- lubridate::ymd(index)
    registry_start_date <- lubridate::ymd(registry_start_date)
    sim_start_date <- index - lubridate::years(max(num_years_to_estimate))

    # NEED SIMULATION:
    #   - have N years > R registry years available
    #   - haven't provided date of death for registry data (not always available)
    if (!((sim_start_date >= registry_start_date) & (!is.null(death_column)))) {

        # Incidence models
        if (is.null(inc_model) & is.null(inc_formula)) {
            stop("Error: Please provide one of inc_model and inc_formula.")
        }
        if (!is.null(inc_model) & !(is.null(inc_formula))) {
            stop("Error: Please provide only one of inc_model and inc_formula.")
        }
        if (is.null(inc_formula )) {
            stop("Error: Functionality for custom incidence objects isn't fully implemented yet. Please provide an 'inc_formula' and use the default homogeneous Poisson process model.")
        }
        if (is.null(inc_model)) {
            inc_model <- fit_exponential_incidence(inc_formula, data)
        }

        # Survival models
        if (!is.null(surv_model) & !(is.null(surv_formula))) {
            stop("Error: Please provide only one of surv_model and surv_formula.")
        }

        if (!is.null(surv_model) & !(missing(dist))) {
            stop("Error: Please provide only one of surv_model and dist.")
        }

        available_dists <- c('lognormal', 'weibull', 'exponential')
        if (!missing(dist) & ! dist %in% available_dists) {
            stop("Error: Please select one of the following distributions: ", paste(available_dists, collapse=', '))
        }

        if (!missing(surv_formula)) {
            surv_model <- build_survreg(surv_formula, data, dist)
        }

        if (missing(surv_model)) {
            stop("Error: Please provide one of surv_model or surv_formula.")
        }

        sim_time <- as.numeric(difftime(index, sim_start_date, units='days'))
        prev_sim <- new_sim_prevalence(data, index, sim_start_date,
                                       sim_time, inc_model, surv_model,
                                       age_col=age_column,
                                       nsims=N_boot)

        # Create column indicating whether contributed to prevalence for each year of interest
        for (year in num_years_to_estimate) {
            # Determine starting incident date
            starting_incident_date <- index - lubridate::years(year)

            # We'll create a new column to hold a binary indicator of whether that observation contributes to prevalence
            col_name <- paste0("prev_", year, "yr")

            # Determine prevalence as incident date is in range and alive at index date
            prev_sim$results[, (col_name) := as.numeric((incident_date > starting_incident_date & incident_date < index) & alive_at_index)]
        }

    } else {
        prev_sim <- NULL
    }

    # Determine point estimates of prevalence by combining simulated and counted values
    names <- sapply(num_years_to_estimate, function(x) paste('y', x, sep=''))
    estimates <- lapply(setNames(num_years_to_estimate, names),
                        new_point_estimate,  # Function
                        prev_sim$results, index, data,
                        counted_formula,
                        registry_start_date,
                        status_column,
                        population_size, proportion, level, precision)

    surv_model <- if (!is.null(prev_sim)) prev_sim$surv_model else NULL
    inc_model <- if (!is.null(prev_sim)) prev_sim$inc_model else NULL

    if (!is.null(counted_formula)) {
        if (!status_column %in% colnames(data)) {
            stop("Error: cannot find status column ", status_column, " in data frame.")
        }
        counted_prev <- new_counted_prevalence(counted_formula, index, data, registry_start_date, status_column)
    } else {
        counted_prev <- NULL
    }
    object <- list(estimates=estimates, simulated=prev_sim$results,
                   counted=counted_prev,
                   surv_model=surv_model,
                   inc_model=inc_model,
                   index_date=index,
                   proportion=proportion,
                   registry_start=registry_start_date,
                   status_col=status_column)

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

new_point_estimate <- function(year, sim_results, index, registry_data, prev_formula, registry_start_date, status_col,
                               population_size=NULL, proportion=100e3,
                               level=0.95, precision=2) {

    # See if need simulation if have less registry data than required
    initial_date <- index - lubridate::years(year)
    need_simulation <- initial_date < registry_start_date

    # Only count prevalence if formula isn't null
    if (!is.null(prev_formula)) {
        count_prev <- new_counted_prevalence(prev_formula, index, registry_data, max(initial_date, registry_start_date), status_col)

        # See if appending prevalence to simulation data or it's entirely counted
        if (initial_date < registry_start_date) {
            stopifnot(!is.null(sim_results))

            col_name <- paste0("prev_", year, "yr")
            sim_contributions <- sim_results[incident_date < registry_start_date][, sum(get(col_name)), by=sim][[2]]  # Return results column
            the_estimate <- count_prev + mean(sim_contributions)

            # Closure to calculate combined standard error
            se_func <- build_se_func(counted_contribs=count_prev, sim_contribs=sim_contributions)

        } else {
            the_estimate <- count_prev
            # Closure to calculate standard error of counted data
            se_func <- build_se_func(counted_contribs=count_prev)
        }
    } else {
        # If don't have counted data then prevalence estimates are entirely simulated
        col_name <- paste0("prev_", year, "yr")
        sim_contributions <- sim_results[, sum(get(col_name)), by=sim][[2]]  # Return results column
        the_estimate <- mean(sim_contributions)

        # Closure to calculate standard error of simulated data
        se_func <- build_se_func(sim_contribs=sim_contributions)
    }

    result <- list(absolute.prevalence=the_estimate)

    if (!is.null(population_size)) {
        the_proportion <- (the_estimate / population_size) * proportion
        se <- se_func(population_size)

        z_level <- qnorm((1+level)/2)
        CI <- z_level * se * proportion

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

build_se_func <- function(counted_contribs=NULL, sim_contribs=NULL) {
    # Pure simulated
    if (is.null(counted_contribs)) {
        function(pop_size) {
            calculate_se_sim(pop_size, sim_contribs)
        }
    }  else if (is.null(sim_contribs)) {
        # Pure counted
        function(pop_size) {
            calculate_se_counted(pop_size, counted_contribs)
        }
    } else {
        # Combination
        function(pop_size) {
            calculate_se_combined(pop_size, counted_contribs, sim_contribs)
        }
    }
}

calculate_se_combined <- function(population_size, counted_contribs, sim_contribs) {
    calculate_se_sim(population_size, sim_contribs) +
        calculate_se_counted(population_size, counted_contribs)
}

calculate_se_sim <- function(population_size, sim_contribs) {
    sd(sim_contribs) / population_size
}
calculate_se_counted <- function(population_size, counted_contribs) {
    raw_proportion <- counted_contribs / population_size
    sqrt((raw_proportion * (1 - raw_proportion)) / population_size)
}

# Start date allows users to specify how long estimating prevalence for, as otherwise including
# all contributions in data set
new_counted_prevalence <- function(formula, index, data, start_date, status_col) {
    death_col <- all.vars(update(formula, .~0))
    entry_col <- all.vars(update(formula, 0~.))

    incident <- data[[entry_col]] >= start_date & data[[entry_col]] < index

    # Use dead at index as it's a simpler boolean operation that just needs negating
    dead_at_index <- !(data[[death_col]] > index) & (data[[status_col]] == 1)

    sum(incident & !dead_at_index)
}

new_sim_prevalence <- function(data, index, starting_date,
                               number_incident_days,
                               inc_model, surv_model,
                               age_col='age',
                               nsims=1000) {

    data <- data[complete.cases(data), ]
    full_data <- data

    covars <- extract_covars(surv_model)

    all_results <- replicate(nsims, {
        # bootstrap dataset
        data <- full_data[sample(seq(nrow(full_data)), replace=T), ]

        # fit incidence and survival models.
        bs_inc <- eval(inc_model$call)
        bs_surv <- eval(surv_model$call)

        # Draw the incident population using the fitted model and predict their death times
        incident_population <- draw_incident_population(bs_inc, full_data, number_incident_days, extract_covars(bs_surv))
        data.table::setDT(incident_population)

        # For each individual determine the time between incidence and the index
        incident_date <- as.Date(starting_date + incident_population$time_to_entry)
        time_to_index <- as.numeric(difftime(index, incident_date, units='days'))

        # Estimate whether alive as Bernouilli trial with p = S(t)
        surv_prob <- predict_survival_probability(bs_surv, incident_population[, -1], time_to_index)
        incident_population[, 'incident_date' := incident_date]
        incident_population[, 'time_to_index' := time_to_index]
        incident_population[, 'alive_at_index' := rbinom(length(surv_prob), size=1, prob=surv_prob)]
        incident_population

    }, simplify=FALSE)

    # Combine into single table
    results <- data.table::rbindlist(all_results, idcol='sim')

    # Force death at 100 if possible
    if (!is.null(age_col) & age_col %in% colnames(results)) {
        results[(get(age_col)*365.25 + time_to_index) > 36525, alive_at_index := 0]
    } else {
        message("No column found for age in ", age_col, ", so cannot assume death at 100 years of age. Be careful of 'infinite' survival times.")
    }

    # This intermediary column isn't useful for the user and would just clutter up the output
    results[, time_to_index := NULL]

    list(results=results,
         surv_model=surv_model,
         inc_model=inc_model)
}

build_curemodel <- function(formula, data, ...) {

}

build_survreg <- function(formula, data, user_dist) {
    # Transforms the registry data into the format specified by survreg.fit,
    # i.e. as a matrix of values with the survival times log transformed.
    complete <- data[complete.cases(data), ]
    X <- model.matrix(formula, complete)
    survobj <- with(complete, eval(formula[[2]]))
    Y <- cbind(survobj[, 1], survobj[, 2])
    data_mat <- cbind(Y, X)

    # Transform of y
    data_mat[, 1] <- survreg.distributions[[user_dist]]$trans(data_mat[, 1])

    # Obtain actual dist
    surv_dist <- survreg.distributions[[user_dist]]$dist
    scale <- if (is.null(survival::survreg.distributions[[user_dist]]$scale)) 0 else survival::survreg.distributions[[user_dist]]$scale

    mod <- survival::survreg.fit(
                         data_mat[, 3:ncol(data_mat)],
                         data_mat[, 1:2],
                         NULL, # weights
                         numeric(nrow(data)), # offset
                         NULL, # init
                         survival::survreg.control(), # controlvars
                         survival::survreg.distributions[[surv_dist]], # dist
                         scale,
                         1, # nstrat
                         0, # strata
                         NULL # parms
                         )

    func_call <- match.call()
    func_call$formula <- eval(formula)
    func_call$user_dist <- eval(user_dist)

    object <- list(coefs=coef(mod),
                   covars = rownames(attr(terms(formula), "factors"))[2:nrow(attr(terms(formula), "factors"))],
                   call = func_call,
                   dist = user_dist,
                   terms= labels(terms(formula))
                   )
    class(object) <- c(class(object), 'survregmin')
    object
}

extract_covars <- function(object) UseMethod("extract_covars")

extract_covars.flexsurvreg <- function(object) {
    attr(object$concat.formula, "covnames")
}

extract_covars.survregmin <- function(object) {
    object$covars
}


predict_survival_probability <- function(object, newdata, times) UseMethod("predict_survival_probability")

predict_survival_probability.survregmin <- function(object, newdata, times) {

    # Expand data into dummy categorical and include intercept
    formula <- as.formula(paste("~", paste(object$terms, collapse='+')))
    wide_df <- model.matrix(formula, newdata)

    # Obtain coefficient for location parameter
    num_params <- distribution_params[[object$dist]]
    if (num_params == 2) {
        lps <- wide_df %*% object$coefs[-length(object$coefs)]
    } else if (num_params == 1) {
        lps <- wide_df %*% object$coefs
    } else {
        stop("Error: Unknown number of parameters ", num_params)
    }

    # Can't see any other way to do this without making subclasses for each distribution
    if (num_params == 2) {
        scale <- object$coefs[length(object$coefs)]
    } else if (num_params == 1) {
        scale <- NULL
    } else {
        stop("Error: Unknown number of parameters ", num_params)
    }

    1- distribution_drawing[[object$dist]](times, lps, scale)
}


# TODO Combine these
distribution_drawing <- list('weibull' = function(times, lps, scale=NULL) {
                                pweibull(times, 1 / exp(scale), exp(lps))
                             },
                             'lognormal' = function(times, lps, scale=NULL) {
                                 plnorm(times, lps, exp(scale))
                             },
                             'loglogistic' = function(times, lps, scale=NULL) {
                                 pllogis(times, 1 / exp(scale), exp(lps))
                             },
                             'exponential' = function(times, lps, scale=NULL) {
                                 pexp(times, 1/exp(lps))
                             }
                             )
distribution_params <- list('weibull'=2,
                            'lognormal'=2,
                            'loglogistic'=2,
                            'exponential'=1
                            )

# TODO Put these in an incidence script somewhere
fit_exponential_incidence <- function(inc_form, data) {
    entry_col <- all.vars(update(inc_form, .~0))
    if (is.null(entry_col) || length(entry_col) != 1) {
        stop("Error: Please provide only a single LHS term in the inc_formula, representing the column holding incident date.")
    }

    registry_duration <- as.numeric(difftime(max(data[[entry_col]]), min(data[[entry_col]]), units='days'))

    strata <- all.vars(update(inc_form, 0~.))
    # Obtain rate for combination of each of these factors and throw error if any instances < required
    if (length(strata) > 0) {

        if (!all(sapply(strata, function(x) is.factor(data[[x]])))) {
            stop("Error: Please format any incidence strata as factors.")
        }

        strata_counts <- table(data[, strata])
        if (any(strata_counts < MIN_INCIDENCE)) {
            stop("Error: Less than ", MIN_INCIDENCE, " incident cases. Please ensure data has sufficient incident cases.")
        }

        rate <- data.frame(strata_counts / registry_duration)
        names(rate)[1:(ncol(rate)-1)] <- strata
        strata_names <- strata
    } else {
        rate <- nrow(data) / registry_duration
        strata_names <- NULL
    }

    # Formulate call with the formula evaluated in this environment with the data argument as a promise
    # that will be evaluated later on the bootstrapped data
    func_call <- match.call()
    func_call$inc_form <- eval(inc_form)

    obj <- list(rate=rate, call=func_call, strata=strata_names)
    class(obj) <- "expinc"
    obj
}

draw_incident_population <- function(object, data, timeframe, covars, ...) UseMethod("draw_incident_population")

draw_incident_population.expinc <- function(object, data, timeframe, covars) {

    if (is.null(object$strata)) {
        initial_num_inds <- INCIDENCE_MARGIN * object$rate * timeframe
        time_to_entry <- cumsum(rexp(initial_num_inds, object$rate))
        time_to_entry <- time_to_entry[time_to_entry < timeframe]
        new_data <- data.table::data.table(time_to_entry)
    } else {
        new_data <- data.table::rbindlist(apply(object$rate, 1, function(row) {
            rate <- as.numeric(row['Freq'])
            initial_num_inds <- INCIDENCE_MARGIN * rate * timeframe
            time_to_entry <- cumsum(rexp(initial_num_inds, rate))
            time_to_entry <- time_to_entry[time_to_entry < timeframe]
            # Form data frame with these factor levels too
            new_data <- data.table::data.table(time_to_entry)
            new_data[, names(row[-length(row)]) := row[-length(row)]]
            new_data
        }))
    }

    # draw covariate values from model that are in survival formula and not in data frame
    missing_covars <- setdiff(covars, colnames(new_data))
    higher_order <- missing_covars[grep("^I\\(.+\\)$", missing_covars)]
    covars_to_sample <- setdiff(missing_covars, higher_order)

    new_covars <- setNames(lapply(covars_to_sample, function(x) sample(data[[x]], nrow(new_data), replace=T)), covars_to_sample)
    data.table::setDT(new_covars)

    # Add in higher order terms
    # NB This method of evaluating higher order won't work if the covariate isn't in this new_covars data, however,
    # this is acceptable as the covariates missing from this data set will be those used in the incidence modelling and
    # so will be categorical and thus not not have higher order terms
    for (term in higher_order) {
        new_covars[, (term) := eval(parse(text=term))]
    }

    # Combine together to form a population that has an entry time and the required covariates for modelling survival
    cbind(new_data, new_covars)
}

predict_survival_probability.flexsurvcure <- function(object, newdata=NULL,
                                                      t=NULL,
                                                      ...)
{
    # Get estimates from flexsurvreg method
    estimates <- predict_survival_probability.flexsurvreg(object, newdata, t, ...)

    if (!is.null(object$pop_mortality)) {
        # Obtain age at index in days
        # TODO Possible to not have hardcoded?
        newdata$age <- floor((newdata$age * 365.25) + t)

        # Obtain mortality rates from these values
        comb <- left_join(newdata, object$pop_mortality, by=c('age', object$pop_covars))

        # For people that don't have a mortality set to 0 as assumedly because they
        # are > 100
        comb[is.na(comb$surv), 'surv'] <- 0

        # scale estimates
        estimates <- estimates * comb$surv
    }

    estimates

}

# This also works for flexsurvcure models
predict_survival_probability.flexsurvreg <- function(object, newdata=NULL,
                                                     t=NULL,
                                                     ...)
{
    x <- object
    dat <- x$data
    Xraw <- model.frame(x)[,unique(attr(model.frame(x),"covnames.orig")),drop=FALSE]
    isfac <- sapply(Xraw, function(x){is.factor(x) || is.character(x)})

    X <- form.model.matrix(object, as.data.frame(newdata))

    fn <- summary.fns(x)
    fn <- expand.summfn.args(fn)  # TODO is this needed?
    beta <- if (x$ncovs==0) 0 else x$res[x$covpars,"est"]

    if (ncol(X) != length(beta)){
        ## don't think we should ever reach here - error should be caught in newdata or X
        isare <- if(length(beta)==1) "is" else "are"
        plural <- if(ncol(X)==1) "" else "s"
        pluralc <- if(length(beta)==1) "" else "s"
        stop("Supplied X has ", ncol(X), " column",plural," but there ",isare," ",
             length(beta), " covariate effect", pluralc)
    }

    dlist <- x$dlist

    # Obtain distribution parameters for each individual
    all_pars <- lapply(1:nrow(X), function(i) {
        basepars.mat <- add.covs(x, x$res.t[dlist$pars,"est"], beta, X[i,,drop=FALSE], transform=FALSE)
        as.list(as.data.frame(basepars.mat))
    })

    # Now convert this to a list of parameters rather than a list of patients
    fnlist <- list(t)
    for (par in dlist$pars) {
        fnlist[[par]] <- sapply(all_pars, function(x) x[[par]])
    }

    do.call(fn, fnlist)
}

form.model.matrix <- function(object, newdata){
    mfo <- model.frame(object)

    ## If required covariate missing, give a slightly more informative error message than, e.g.
    ## "Error in eval(expr, envir, enclos) (from flexsurvreg.R#649) : object 'sex' not found"
    covnames <- attr(mfo, "covnames")
    missing.covs <- unique(covnames[!covnames %in% names(newdata)])
    if (length(missing.covs) > 0){
        missing.covs <- sprintf("\"%s\"", missing.covs)
        plural <- if (length(missing.covs)>1) "s" else ""
        stop(sprintf("Value%s of covariate%s ",plural,plural), paste(missing.covs, collapse=", "), " not supplied in \"newdata\"")
    }

    ## as in predict.lm
    tt <- attr(mfo, "terms")
    Terms <- delete.response(tt)
    mf <- model.frame(Terms, newdata, xlev = .getXlevels(tt, mfo))
    if (!is.null(cl <- attr(Terms, "dataClasses")))
        .checkMFClasses(cl, mf)

    forms <- object$all.formulae
    mml <- vector(mode="list", length=length(object$dlist$pars))
    names(mml) <- names(forms)
    forms[[1]] <- delete.response(terms(forms[[1]]))
    for (i in names(forms)){
        mml[[i]] <- model.matrix(forms[[i]], mf)
    }
    X <- compress.model.matrices(mml)

    attr(X, "newdata") <- mf # newdata with any extra variables stripped.  Used to name components of summary list
    X
}

compress.model.matrices <- function(mml){
    cbind.drop.intercept <- function(...)do.call("cbind", lapply(list(...), function(x)x[,-1,drop=FALSE]))
    X <- do.call("cbind.drop.intercept",mml)
    loc.cnames <- colnames(mml[[1]])[-1]
    anc.cnames <- unlist(mapply(function(x,y)sprintf("%s(%s)",x,y), names(mml[-1]), lapply(mml[-1], function(x)colnames(x)[-1])))
    cnames <- c(loc.cnames, anc.cnames)
    colnames(X) <- cnames
    X
}

summary.fns <- function(x){
   function(t, ...) {
       1 - x$dfns$p(t,...)
   }
}

expand.summfn.args <- function(summfn){
    summfn2 <- summfn
    args <- c(alist(t=), formals(summfn))
    formals(summfn2) <- args[!duplicated(names(args))]
    body(summfn2) <- body(summfn)
    summfn2
}

add.covs <- function(x, pars, beta, X, transform=FALSE){  ## TODO option to transform on input
    nres <- nrow(X)
    if (!is.matrix(pars)) pars <- matrix(pars, nrow=nres, ncol=length(pars), byrow=TRUE)
    if (!is.matrix(beta)) beta <- matrix(beta, nrow=1)
    for (j in seq(along=x$dlist$pars)){
        covinds <- x$mx[[x$dlist$pars[j]]]
        if (length(covinds) > 0){
            pars[,j] <- pars[,j] + beta[,covinds] %*% t(X[,covinds,drop=FALSE])
        }
        if (!transform)
            pars[,j] <- x$dlist$inv.transforms[[j]](pars[,j])
    }
    colnames(pars) <- x$dlist$pars
    pars
}
