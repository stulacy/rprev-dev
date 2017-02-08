# rprev 0.2.3
Incidence function returning different values to the known_inc_rate function from a prevalence object
This is due to how we define years. The prevalence function calculates a year based on 365.25 days and calculates the starting date of the time-period of interest using the above method. Whereas in your manual call to incidence you are interested in incidence from the 9 years from 2005-09-01 to 2014-09-01, and so specify the starting date as 2005-09-01. However, due to the extra .25 of a day, 9 years back from 2014-09-01 is actually 2005-08-31, so any patient incident on the first of September in any year gets shifted years. I'm happy to keep using our definition of years, although in order to aid user friendliness I've reparameterised incidence (see below) to allow either the starting date, the ending date, or both, similarly to the seq function in base R.

Renamed raw_incidence to yearly_incidence and reparameterised
Renamed to be more descriptive of what the function actually does, and reparameterised to allow the user to specify the ending date of the time interval of interested instead. Note I haven't removed raw_incidence, just deprecated it in favour of yearly_incidence.

Renamed determine_registry_years to determine_yearly_limits and reparameterised
I realised that the original function name isn't very descriptive for what it does (provides the yearly end points of a specific time interval) and so have tried to reflect its purpose with the new name. Not entirely sold on it so feel free to suggest improvements! Also have reparameterised so the user can specify the closing date in the interval rather than the opening.

 plot.incidence not returning ggplot object
It now does (along with other plot methods), allowing for easier manual tweaking. I'm not sure why I didn't have this behaviour originally, there must have been a reason...

calculating prevalence rate / 100K and CIs from counted data
This functionality was always present, but I've now made it more explicit by **not running the simulation when there is more registry data available than needed to estimate N-year prevalence** I'm a bit worried by this as why didn't I think of this before? The previous warning messaged remarked that survival models were still built using all available registry data, but these aren't required to count prevalence. I won't push this change through to release until I've checked with Simon that I'm not completely breaking some theoretical grounds, but in terms of implementation it actually changes nothing.

ONS citation:
Updated to include a link to the specific webpage where the data set is obtained from and improved formatting.

other:

    summary.prevalence correctly displaying posterior age distributions of simulated cases
    summary.prevalence displaying the prevalence estimates themselves
    unit tests updated to reflect the above changes
    vignette updated to reflect the above changes



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
