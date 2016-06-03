#' rprev: Estimate disease point prevalence using a combination of registry data and
#' Monte Carlo simulations.
#'
#' The rprev package uses available registry data to estimate point prevalence
#' at a specified index date. This is done by calculating yearly incident cases
#' and estimating survival probabilities of these cases at the index date, to
#' establish yearly contributions of incidence to the point prevalence estimate.
#'
#' This process relies upon accurate modeling of both the incidence and survival
#' process, requiring that two assumptions are met: \itemize{ \item That the
#' disease incidence is a homogeneous Poisson process \item That survival can be
#' modeled by a Weibull model, incorporating both age and sex }
#'
#' Prevalence is estimated using incident cases from a set number of years,
#' where the larger this values the more accurate the prevalence estimates are.
#' However, if the user asks to use more years of incident cases than are
#' available in the registry data set, then the remaining years of incidence are
#' simulated.
#'
#' The primary function in this package is thereby \code{\link{prevalence}},
#' which performs the combination of counted incidence from the registry data,
#' and the simulated cases, along with the calculation of their survival
#' probabilities at the index date.
#'
#' \code{\link{incidence}} provides a summary of the incident cases in the
#' registry data set, allowing for inspection of whether the homogeneous Poisson
#' process assumption holds for the disease in question.
#'
#' @docType package
#' @name rprev
NULL
