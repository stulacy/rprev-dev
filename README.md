# rprev

rprev estimates disease prevalence at a specified index date from registry data. To improve the estimate accuracy, Monte Carlo simulation techniques are used to simulate incident cases in years for which incidence data is unavailable. Prevalence arises from two independent stochastic processes: disease incidence and survival. This package provides functionality for providing models for these processes in an object-oriented manner. Default models are provided that model incidence as a homogeneous Poisson process and survival as a standard parameteric distribution. See the *user_guide* vignette for more details about the implementation, and the original publication for details of the algorithm, available at [http://www.ncbi.nlm.nih.gov/pubmed/24656754](http://www.ncbi.nlm.nih.gov/pubmed/24656754).

## Installation

To install from CRAN, simply use `install.packages('rprev')`, while the latest development version can be installed from GitHub using `devtools::install_github("stulacy/rprev-dev")`.
