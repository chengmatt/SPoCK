.onAttach <- function(libname, pkgname) {
  # Check for OSA_multivariate_dists
  if (!requireNamespace("OSA_multivariate_dists", quietly = TRUE)) {
    message("Installing OSA_multivariate_dists...")
    tryCatch({
      TMB:::install.contrib("https://github.com/vtrijoulet/OSA_multivariate_dists/archive/main.zip")
    }, error = function(e) {
      warning("Failed to install OSA_multivariate_dists: ", conditionMessage(e))
    })
  }

  # Check for compResidual
  if (!requireNamespace("compResidual", quietly = TRUE)) {
    message("Installing compResidual from GitHub...")
    tryCatch({
      remotes::install_github("fishfollower/compResidual/compResidual", force = TRUE)
    }, error = function(e) {
      warning("Failed to install compResidual: ", conditionMessage(e))
    })
  }
}
