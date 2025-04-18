# install OSA_multivariate_dists
install_osa_multivariate <- function() {
  if (!requireNamespace("TMB", quietly = TRUE)) {
    stop("TMB package required but not installed.")
  }

  # call unexported TMB function
  get("install.contrib", envir = asNamespace("TMB"))(
    "https://github.com/vtrijoulet/OSA_multivariate_dists/archive/main.zip"
  )
}

# install compResidual
install_compResidual <- function() {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    stop("remotes package required but not installed.")
  }

  tryCatch({
    remotes::install_github("fishfollower/compResidual/compResidual", force = TRUE)
  }, error = function(e) {
    stop("Failed to install compResidual: ", conditionMessage(e))
  })
}

.onAttach <- function(libname, pkgname) {
  # check and install OSA_multivariate_dists if not already installed
  if (!requireNamespace("OSA_multivariate_dists", quietly = TRUE)) {
    message("Installing OSA_multivariate_dists from GitHub...")
    tryCatch({
      install_osa_multivariate()
    }, error = function(e) {
      packageStartupMessage("Failed to install OSA_multivariate_dists: ", conditionMessage(e))
    })
  }

  # check and install compResidual if not already installed
  if (!requireNamespace("compResidual", quietly = TRUE)) {
    message("Installing compResidual from GitHub...")
    tryCatch({
      install_compResidual()
    }, error = function(e) {
      packageStartupMessage("Failed to install compResidual: ", conditionMessage(e))
    })
  }
}
