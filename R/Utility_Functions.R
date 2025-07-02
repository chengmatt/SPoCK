#' ggplot theme for sablefish
#'
#' @return ggplot theme
#' @export theme_sablefish
#' @import ggplot2
theme_sablefish <- function() {
   theme_bw() +
    theme(legend.position = "top",
          strip.text = element_text(size = 17),
          title = element_text(size = 21, color = 'black'),
          axis.text = element_text(size = 15, color = "black"),
          axis.title = element_text(size = 17, color = 'black'),
          legend.text = element_text(size = 15, color = "black"),
          legend.title = element_text(size = 17, color = 'black'))
}


#' Function to fill in an n x n correlation AR(1) matrix
#'
#' @param n Number of bins
#' @param rho correaltion parameter
#'
#' @return correlation matrix for an ar1 process
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' mat <- get_AR1_CorrMat(10, 0.5)
#' }
get_AR1_CorrMat <- function(n, rho) {
  corrMatrix <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      # Calculate the correlation based on the lag distance
      corrMatrix[i, j] <- rho^(abs(i - j))
    } # end i
  } # end j
  return(corrMatrix)
}

#' Constant correlation matrix
#'
#' @param n Number of bins
#' @param rho correaltion parameter
#'
#' @return constant correlation matrix
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' mat <- get_Constant_CorrMat(10, 0.5)
#' }
get_Constant_CorrMat <- function(n, rho) {
  corrMatrix <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      if(i != j) corrMatrix[i, j] <- rho
      else corrMatrix[i, j] <- 1
    } # end i
  } # end j
  return(corrMatrix)
}

#' For combining a parameter and data list in RTMB so a data object can be explicitly defined
#'
#' @param f Parameter list
#' @param d Data list
#' @keywords internal
#' @examples
#' \dontrun{
#'   obj <- RTMB::MakeADFun(cmb(sabie_RTMB, data), parameters = parameters, map = mapping, random = random, silent = TRUE)
#' }
cmb <- function(f, d) {
  function(p) f(p, d)
}

#' Helper function to collect messages
#'
#' @param ... character vector of messages
#' @keywords internal
collect_message <- function(...) {
  messages_list <<- c(messages_list, paste(..., sep = ""))
}

#' Helper function to check package availbility
#'
#' @param pkg package name character
#' @keywords internal
is_package_available <- function(pkg) {
  nzchar(system.file(package = pkg))
}

#' Go from TAC to Fishing Mortality using bisection
#'
#' @param f_guess Initial guess of F
#' @param catch Provided catch values
#' @param NAA Numbers, dimensioned by ages, and sexes
#' @param WAA Weight, dimensioned by ages and sexes
#' @param natmort Natural mortality dimensioned by ages and sex
#' @param fish_sel Fishery selectivity, dimesnioned by ages and sex
#' @param n.iter Number of iterations for bisection
#' @param lb Lower bound of F
#' @param ub Upper bound of F
#'
#' @returns Fishing mortality values
#' @export bisection_F
bisection_F <- function(f_guess,
                        catch,
                        NAA,
                        WAA,
                        natmort,
                        fish_sel,
                        n.iter = 20,
                        lb = 0,
                        ub = 2) {

  range <- vector(length=2) # F range
  range[1] <- lb # Lower bound
  range[2] <- ub # Upper bound

  for(i in 1:n.iter) {

    # Get midpoint of range
    midpoint <- mean(range)

    # Caclulate baranov's
    FAA <- (midpoint * fish_sel)
    ZAA <- FAA + natmort
    pred_catch <- sum((FAA / ZAA * NAA * (1 - exp(-ZAA))) * WAA)

    if(pred_catch < catch) {
      range[1] <- midpoint
      range[2] <- range[2]
    }else {
      range[1] <- range[1]
      range[2] <- midpoint
    }

  } # end i loop

  return(midpoint)
}

#' Post Optimization Model Convergence Checks
#'
#' @param sd_rep sd report list from a `SPoCK` model
#' @param rep report list from a `SPoCK` model
#' @param gradient_tol Value for maximum gradient tolerance to use
#' @param se_tol Value for maximum standard error tolerance to use
#' @param corr_tol Value for maximum correlation tolerance to use
#'
#' @export post_optim_sanity_checks
post_optim_sanity_checks <- function(sd_rep,
                                     rep,
                                     gradient_tol = 1e-3,
                                     se_tol = 100,
                                     corr_tol = 0.99
                                     ) {

  passed_post_sanity_checks <- TRUE

  # check likelihoods are all finite and not NA
  if(!all(is.finite(rep$jnLL))) {
    message("Found Inf in joint log-likelihood, model is not converged!")
    passed_post_sanity_checks <- F
  }

  # check maximum absolute gradients
  max_abs_grad_ndx <- which.max(abs(sd_rep$gradient.fixed))
  max_abs_grad <- abs(sd_rep$gradient.fixed)[max_abs_grad_ndx]
  if(gradient_tol < max_abs_grad) {
    message("Parameter: ", names(sd_rep$par.fixed)[max_abs_grad_ndx], " had absolute gradient = ", max_abs_grad,
            " which was greater than tolerance ", gradient_tol,". This indicates potential non-convergence according to the tolerance.\n")
    passed_post_sanity_checks <- F
  }

  # check hessian
  if(!sd_rep$pdHess) {
    message("Hessian is not positive definite, model is not converged!")
    passed_post_sanity_checks <- F
  }

  # check if standard errors are finite
  if(!all(is.finite(sqrt(diag(sd_rep$cov.fixed))))) {
    message("Found non finite elements in standard errors of parameters, model is not converged!")
    passed_post_sanity_checks <- F
  }

  # check if standard errors are big
  if(max(diag(sd_rep$cov.fixed)) > se_tol) {
    message("Parameter: ", names(diag(sd_rep$cov.fixed))[which.max(diag(sd_rep$cov.fixed))], " has a standard error = ",
            max(diag(sd_rep$cov.fixed)), " which was greated than tolerance ", se_tol, ". This indicates potential non-convergence according to the tolerance. \n")
    passed_post_sanity_checks <- F
  }

  # check if correlations are big
  corr_mat <- cov2cor(sd_rep$cov.fixed)
  diag(corr_mat) <- "Same" # set diagonal to "Same" to remove from max calculations

  # reshape to dataframe
  corr_df <- reshape2::melt(corr_mat) %>%
    dplyr::filter(value != 'Same') %>%
    dplyr::mutate(value = as.numeric(value))

  if(max(abs(corr_df$value)) > corr_tol) {
    message("Parameter pairs: ", corr_df$Var1[which.max(abs(corr_df$value))], " and ", corr_df$Var2[which.max(abs(corr_df$value))], " have a correlation of ", max(abs(corr_df$value)), ". This indicates potential non-convergence according to the tolerance.")
    passed_post_sanity_checks <- F
  }

  cat("\n\n");
  if(passed_post_sanity_checks) {
    message("Successfully passed post-optim-sanity checks\n")
    return(TRUE)
  }

}
