# rprev 0.2.0

Minor bug fixes and a slight change to the parameterisation of prevalence:

  - In prevalence, prevalence_counted, and prevalence_simulated, the user specifies the index date at which to estimate prevalence, rather than having it inferred from the data
  - max_yearly_incidence has been removed as a parameter from both prevalence and prevalence_simulated as it can be calculated from the supplied data
  
# rprev 0.1.0

First release of the package, working with all features necessary to provide estimates of point prevalence. Issues which we'd like to address in future releases are:

  - Allow for other incidence processes than homogeneous Poisson
  - Enable more flexibility in survival modelling, rather than Weibull regression with linear covariate effects
  - Allow for the inclusion of more covariates in both the survival modelling, and the marking of the incidence process
