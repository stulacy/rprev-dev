#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}

// [[Rcpp::export]]
NumericMatrix bootstrapfit(DataFrame) {

}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
timesTwo(42)
  form = Surv(survival_time, indicator) ~ age_initial + sex
  data = data_r
  n = N_boot

  orig = survreg2(form, data=data)
  coef(orig)
  orig$scale

  complete = data[complete.cases(data),]
  X = model.matrix(form, data)

  # TRANSFORM Y
  Y = as.matrix(complete[, c('survival_time', 'indicator')])
	tranfun <- survreg.distributions[['weibull']]$trans
	Y2 <- cbind(tranfun(Y[,1]), Y[,2])

  test = survreg.fit(X,
                     Y2,
                     NULL, # weights
                     numeric(nrow(complete)), # offset
                     NULL, # init
                     survreg.control(), # controlvars
                     survreg.distributions[['extreme']], # dist
                     0, # scale
                     1, # nstrat
                     0, # strata
                     NULL # parms
                     )
  coef(test)[1:ncol(X)]
  exp(coef(test)[-(1:ncol(X))])

              #ith_data <- data[sample(1:(dim(data)[1]), dim(data)[1], replace=T), ]
              #wbb <- survreg(form, data=ith_data)
              #c(wbb$coe, wbb$scale)
*/
