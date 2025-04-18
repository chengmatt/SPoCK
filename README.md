# SPoCK: Spatial (S) Population (Po) Model with Connectivity (C) across generalized Komponents (K)
A generalized spatial stock assessment model written in RTMB, than be extended to represent any number of regions, ages, sexes, and fleets. 

## Installation

`SPoCK` is a package written in RTMB, and relies on several packages for plotting and model diagnostics purposes. To install the package, users should have `devtools` installed. It is also generally recommended to install the packaged listed below in the following order

```
install.packages("devtools")
install.packages("TMB")
install.packages("RTMB")
TMB:::install.contrib("https://github.com/vtrijoulet/OSA_multivariate_dists/archive/main.zip")
devtools::install_github("fishfollower/compResidual/compResidual")
devtools::install_github("noaa-afsc/afscOSA", dependencies = TRUE)
devtools::install_github("chengmatt/SPoCK", dependencies = TRUE)
```
Some users may experience installation issues regarding permissions. This can potentially be circumvented with the following code:
```
devtools::install_github("chengmatt/SPoCK", dependencies = TRUE, lib = Sys.getenv("R_LIBS_USER"))
```

