build_survreg <- function(formula, data, user_dist) {
    # Transforms the registry data into the format specified by survreg.fit,
    # i.e. as a matrix of values with the survival times log transformed.
    complete <- data[complete.cases(data), ]
    X <- model.matrix(formula, complete)
    survobj <- with(complete, eval(formula[[2]]))
    Y <- cbind(survobj[, 1], survobj[, 2])

    # Obtain covariate names
    factors <- attr(terms(formula), "factors")
    covars <- if (length(factors) == 0) {
                    c()
               } else {
                   row.names(factors)[2:nrow(factors)]
               }

    # Transform of y
    dist_obj <- survival::survreg.distributions[[user_dist]]
    Y[, 1] <- dist_obj$trans(Y[, 1])

    # Obtain actual dist
    surv_dist <- dist_obj$dist
    scale <- if (is.null(dist_obj$scale)) 0 else dist_obj$scale


    mod <- survival::survreg.fit(
                         X,
                         Y,
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
                   covars = covars,
                   call = func_call,
                   dist = user_dist,
                   terms= labels(terms(formula))
                   )
    class(object) <- c(class(object), 'survregmin')
    object
}

#' @export
extract_covars.survregmin <- function(object) {
    object$covars
}

#' @export
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
