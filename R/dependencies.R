# Internal variables to track dependency availability
.pkg_env <- new.env(parent = emptyenv())

#' Install SPoCK dependencies
#'
#' This function installs the optional dependencies required for full
#' functionality of the SPoCK package.
#'
#' @param force Logical, whether to force reinstallation of dependencies even if they are already installed
install_spock_dependencies <- function(force = FALSE) {
  # Install OSA_multivariate_dists (TMB extension module)
  install_osa_multivariate <- function() {
    message("Installing OSA_multivariate_dists TMB extension...")
    if (!requireNamespace("TMB", quietly = TRUE)) {
      stop("TMB package required but not installed.")
    }
    tryCatch({
      get("install.contrib", envir = asNamespace("RTMB"))(
        "https://github.com/vtrijoulet/OSA_multivariate_dists/archive/main.zip"
      )
      message("OSA_multivariate_dists TMB extension successfully installed.")
      TRUE
    }, error = function(e) {
      warning("Failed to install OSA_multivariate_dists: ", conditionMessage(e))
      FALSE
    })
  }

  # For OSA module, we need a different check since it's not a package
  osa_installed <- tryCatch({
    # Check if we can find the TMB module
    tmb_path <- file.path(system.file("include", package = "TMB"), "contrib")
    any(dir.exists(file.path(tmb_path, "OSA_multivariate_dists-main")))
  }, error = function(e) FALSE)

  if (force || !osa_installed) {
    .pkg_env$has_osa <- install_osa_multivariate()
  } else {
    .pkg_env$has_osa <- TRUE
  }

  # Install compResidual
  if (force || !requireNamespace("compResidual", quietly = TRUE)) {
    message("Installing compResidual from GitHub...")
    if (!requireNamespace("remotes", quietly = TRUE)) {
      stop("remotes package required but not installed.")
    }
    tryCatch({
      remotes::install_github("fishfollower/compResidual/compResidual", force = TRUE)
    }, error = function(e) {
      warning("Failed to install compResidual: ", conditionMessage(e))
    })
  }

  # Install afscOSA
  if (force || !requireNamespace("afscOSA", quietly = TRUE)) {
    message("Installing afscOSA from GitHub...")
    if (!requireNamespace("remotes", quietly = TRUE)) {
      stop("remotes package required but not installed.")
    }
    tryCatch({
      remotes::install_github("noaa-afsc/afscOSA", force = TRUE)
    }, error = function(e) {
      warning("Failed to install afscOSA: ", conditionMessage(e))
    })
  }

  # Update availability status
  check_dependencies()

  # Return availability status
  return(list(
    OSA_multivariate_dists = .pkg_env$has_osa,
    compResidual = .pkg_env$has_compResidual,
    afscOSA = .pkg_env$has_afscOSA
  ))
}

# Internal function to check dependency availability
check_dependencies <- function() {
  # For OSA module, check if the directory exists
  .pkg_env$has_osa <- tryCatch({
    tmb_path <- file.path(system.file("include", package = "TMB"), "contrib")
    any(dir.exists(file.path(tmb_path, "OSA_multivariate_dists-main")))
  }, error = function(e) FALSE)

  .pkg_env$has_compResidual <- requireNamespace("compResidual", quietly = TRUE)
  .pkg_env$has_afscOSA <- requireNamespace("afscOSA", quietly = TRUE)
}

# Replace the problematic .onAttach with .onLoad
.onLoad <- function(libname, pkgname) {
  # Check which dependencies are available
  check_dependencies()
}

.onAttach <- function(libname, pkgname) {
  # Inform users about dependency status
  if (!.pkg_env$has_osa) {
    packageStartupMessage("OSA_multivariate_dists TMB extension is not installed. Some analytical functions will be limited (OSA stuff).")
    packageStartupMessage("Run install_spock_dependencies() to install all required dependencies.")
  }

  if (!.pkg_env$has_compResidual) {
    packageStartupMessage("compResidual is not installed. Model diagnostics will be limited (OSA stuff).")
  }

  if (!.pkg_env$has_afscOSA) {
    packageStartupMessage("afscOSA is not installed. Some AFSC-specific functionality will be unavailable (OSA stuff).")
  }
}

#' Check availability of dependency packages
#'
#' This function returns the availability status of all optional dependencies
#' used by the SPoCK package.
#'
#' @return A named list with logical values indicating if each dependency is available
has_dependencies <- function() {
  check_dependencies() # Refresh status
  return(list(
    OSA_multivariate_dists = .pkg_env$has_osa,
    compResidual = .pkg_env$has_compResidual,
    afscOSA = .pkg_env$has_afscOSA
  ))
}
