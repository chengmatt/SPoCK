# Install optional dependencies, but leave it to the user to decide
install_spock_dependencies <- function() {
  # Check if OSA_multivariate_dists is installed
  osa_installed <- tryCatch({
    # Check if we can find the TMB module
    tmb_path <- file.path(system.file("include", package = "TMB"), "contrib")
    any(dir.exists(file.path(tmb_path, "OSA_multivariate_dists-main")))
  }, error = function(e) FALSE)

  # Check if compResidual is installed
  compResidual_installed <- requireNamespace("compResidual", quietly = TRUE)

  # Check if afscOSA is installed
  afscOSA_installed <- requireNamespace("afscOSA", quietly = TRUE)

  # Notify the user about missing dependencies
  if (!osa_installed) {
    message("OSA_multivariate_dists TMB extension is not installed. You can install it by running:\n")
    message("  install.packages('TMB')")
    message("  remotes::install_github('vtrijoulet/OSA_multivariate_dists')")
  }

  if (!compResidual_installed) {
    message("compResidual is not installed. You can install it by running:\n")
    message("  remotes::install_github('fishfollower/compResidual/compResidual')")
  }

  if (!afscOSA_installed) {
    message("afscOSA is not installed. You can install it by running:\n")
    message("  remotes::install_github('noaa-afsc/afscOSA')")
  }

  # Return the installation status
  return(list(
    OSA_multivariate_dists = osa_installed,
    compResidual = compResidual_installed,
    afscOSA = afscOSA_installed
  ))
}

.onAttach <- function(libname, pkgname) {
  # Check for dependency availability
  deps <- install_spock_dependencies()

  # Inform users about missing dependencies
  if (!deps$OSA_multivariate_dists) {
    packageStartupMessage("OSA_multivariate_dists TMB extension is not installed. Analytical functions may be limited.")
  }

  if (!deps$compResidual) {
    packageStartupMessage("compResidual is not installed. Some model diagnostics will be unavailable.")
  }

  if (!deps$afscOSA) {
    packageStartupMessage("afscOSA is not installed. AFSC-specific functionality may be unavailable.")
  }
}

