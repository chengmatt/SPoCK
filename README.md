# SPoCK: Spatial (S) Population (Po) Model with Connectivity (C) across generalized Komponents (K)
#### Matt Cheng, Dan Goethel, Pete Hulson, Curry Cunningham

A generalized stock assessment model written in RTMB, that can be extended to represent any number of regions, ages, sexes, and fleets. 

## Installation

`SPoCK` is a package written in RTMB and optionally relies on several packages for plotting and model diagnostics purposes. To install the package, users should have `devtools` installed. It is also generally recommended to install the packages listed below, but need not be if users do not want functionality with one-step ahead residuals.

```
install.packages("devtools") # install dev tools
install.packages("TMB") # install TMB
install.packages("RTMB") # install RTMB
TMB:::install.contrib("https://github.com/vtrijoulet/OSA_multivariate_dists/archive/main.zip") # get multivariate OSA distributions

# optional packages to install
devtools::install_github("fishfollower/compResidual/compResidual") 
devtools::install_github("noaa-afsc/afscOSA", dependencies = TRUE)

# install SPoCK
devtools::install_github("chengmatt/SPoCK", dependencies = TRUE)
```

Some users may experience installation issues regarding permissions. This can potentially be circumvented with the following code:
```
devtools::install_github("chengmatt/SPoCK", dependencies = TRUE, lib = Sys.getenv("R_LIBS_USER"))
```
_This package is still very much under active development. Use at your own risk!_

