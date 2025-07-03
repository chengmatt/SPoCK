# SPoRC: A Generalized Stochastic (S) Population (Po) Model over Regional (R) Components (C)
#### Matt Cheng, Dan Goethel, Pete Hulson, Curry Cunningham

`SPoRC` is a flexible modeling framework that captures population dynamics across space, incorporating stochasticity in vital rates and movement among geographically defined components. It supports integration of multiple data sources, regional structuring, and age- and sex-specific processes, making it well suited for complex metapopulation or spatial stock assessment contexts.

## Installation

`SPoRC` is a package written in RTMB and optionally relies on several packages for plotting and model diagnostics purposes. To install the package, users should have `devtools` installed. It is also generally recommended to install the packages listed below, but need not be if users do not want functionality with one-step ahead residuals.

```
install.packages("devtools") # install dev tools
install.packages("TMB") # install TMB
install.packages("RTMB") # install RTMB
TMB:::install.contrib("https://github.com/vtrijoulet/OSA_multivariate_dists/archive/main.zip") # get multivariate OSA distributions

# optional packages to install
devtools::install_github("fishfollower/compResidual/compResidual") 
devtools::install_github("noaa-afsc/afscOSA", dependencies = TRUE)

# install SPoRC
devtools::install_github("chengmatt/SPoRC", dependencies = TRUE)
```

Some users may experience installation issues regarding permissions. This can potentially be circumvented with the following code:
```
devtools::install_github("chengmatt/SPoRC", dependencies = TRUE, lib = Sys.getenv("R_LIBS_USER"))
```
_This package is still very much under active development. Use at your own risk!_

