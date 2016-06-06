# rprev

rprev estimates disease prevalence at a specified index date from registry data. To improve the estimate accuracy, Monte Carlo simulation techniques are used to simulate incident cases in years for which incidence data is unavailable. Disease survival is modelled with parametric Weibull regression. See the *user_guide* vignette for more details about the implementation, and the original publication for details of the algorithm, available at [http://www.ncbi.nlm.nih.gov/pubmed/24656754](http://www.ncbi.nlm.nih.gov/pubmed/24656754).

## Installation

To install from CRAN, simply use `install.packages('rprev')`. The code is currently not available elsewhere.
