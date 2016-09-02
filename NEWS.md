# rprev 0.2.2

Bug hotfix.

# rprev 0.2.1

The posterior age distribution, returned from `prevalence` as in the `simulated` object, is now stored in the format of a nested list rather than a matrix as before. The first dimension of the list corresponds to each sex (if applicable), the next indexing the number of years of simulated cases, and the final corresponds to the bootstrap samples. The final level comprises a vector holding the ages of the simulated cases which are still contributing to prevalence at the index date from the corresponding sex, year, and bootstrap sample number.

# rprev 0.2.0

Minor bug fixes and a slight change to the parameterisation of prevalence:

  - In prevalence, prevalence_counted, and prevalence_simulated, the user specifies the index date at which to estimate prevalence, rather than having it inferred from the data
  - max_yearly_incidence has been removed as a parameter from both prevalence and prevalence_simulated as it can be calculated from the supplied data
  - prevalence per 100K estimates now have the confidence intervals the correct way around
  - unit tests for prevalence functions don't rely on cached results any longer. This has helped to reduce the size of the source code from 25MB to 2MB.
  
# rprev 0.1.0

First release of the package, working with all features necessary to provide estimates of point prevalence. Issues which we'd like to address in future releases are:

  - Allow for other incidence processes than homogeneous Poisson
  - Enable more flexibility in survival modelling, rather than Weibull regression with linear covariate effects
  - Allow for the inclusion of more covariates in both the survival modelling, and the marking of the incidence process
