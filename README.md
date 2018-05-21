# rprev

rprev estimates disease prevalence at a specified index date from incomplete registry data.
It is designed to be used when estimates of point prevalence from registry data are required, but the registry hasn't been running for sufficiently long to count the number of prevalent cases.
Monte Carlo simulation techniques are used to simulate incident cases in years for which incidence data is unavailable, and then estimate survival at the index date.

Prevalence arises from two independent stochastic processes: disease incidence and survival.
Default models are provided that model incidence as a homogeneous Poisson process and survival as a standard parameteric distribution, although both of these models can be user specified for further control.
See the *user_guide* vignette for more details about the implementation, and the original publication for details of the algorithm, available at [http://www.ncbi.nlm.nih.gov/pubmed/24656754](http://www.ncbi.nlm.nih.gov/pubmed/24656754).

## Installation

To install from CRAN, simply use `install.packages('rprev')`, while the latest development version can be installed from GitHub using `devtools::install_github("stulacy/rprev-dev")`.
