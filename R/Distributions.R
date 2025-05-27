#' Symmetric Beta Function
#'
#' @param p_val Parameter value
#' @param p_ub Upper Bound of Parameter
#' @param p_lb Lower Bound of Parameter
#' @param p_prsd SD of parameter, higher values have a stronger penalty on bounds, lower values have a more difuse penalty on bounds
#' @param log whether or not to return the log likelihood
#'
#' @returns Returns likelihood values from a symmetric beta
#' @keywords internal
#'
dbeta_symmetric <- function(p_val, p_ub, p_lb, p_prsd, log = TRUE) {
  # Calculate mu term
  mu <- p_prsd * log((p_ub + p_lb)/2 - p_lb) - p_prsd * log(0.5)
  # Calculate the prior likelihood components
  term1 <- mu
  term2 <- p_prsd * log(p_val - p_lb + 1e-4)
  term3 <- p_prsd * log(1 - (p_val - p_lb - 1e-4)/(p_ub - p_lb))
  # Combine terms to get final prior likelihood
  nLL <- term1 + term2 + term3
  if(log == TRUE) return(nLL) else return(exp(nLL))
}


#' Dirichlet Likelihood
#'
#' @param x vector of values
#' @param alpha Expected values w/ concentration sum(alpha)
#' @param log Whether to give log or not
#'
#' @returns Returns likelihood values form a dirichlet
#' @keywords internal
#'
ddirichlet <- function(x, alpha, log = TRUE) {
  logres = lgamma(sum(alpha)) - sum(lgamma(alpha)) + sum((alpha - 1) * log(x))
  if(log == TRUE) return(logres) else return(exp(logres))
} # end function

#' Dirichlet Mutlinomial Likelihood
#' From https://github.com/James-Thorson/CCSRA/blob/main/inst/executables/CCSRA_v9.cpp
#' @param obs Vector of observed values in proportions
#' @param pred Vector or predicted values in proportions
#' @param Ntotal Input sample size scalar
#' @param ln_theta Weighting parameter in log space
#' @param give_log Whether or not likelihood is in log space
#'
#' @returns returns likelihood values from a dirihclet multinomial
#' @keywords internal
#'
ddirmult = function(obs, pred, Ntotal, ln_theta, give_log = TRUE) {
  # Set up function variables
  n_c = length(obs) # number of categories
  p_exp = pred # expected values container
  p_obs = obs # observed values container
  dirichlet_Parm = exp(ln_theta) * Ntotal # Dirichlet alpha parameters

  # set up pdf
  logres = lgamma(Ntotal + 1)
  for(c in 1:n_c) logres = logres - lgamma(Ntotal*p_obs[c]+1) # integration constant
  logres = logres + lgamma(dirichlet_Parm) - lgamma(Ntotal+dirichlet_Parm) # 2nd term in formula

  # Summation in 3rd term in formula
  for(c in 1:n_c) {
    logres = logres + lgamma(Ntotal*p_obs[c] + dirichlet_Parm*p_exp[c])
    logres = logres - lgamma(dirichlet_Parm * p_exp[c])
  } # end c

  if(give_log == TRUE) return(logres)
  else return(exp(logres))
} # end function

#' Logistic Normal Likelihood
#'
#' @param obs Vector of observed values in numbers (can be integers or non-integers - vector length of n_bins)
#' @param pred Vector of predicted values in proportions (vector length of n_bins)
#' @param Sigma Covariance matrix
#' @param give_log whether or not to use log space likelihood
#' @param perc percentage to multiply minimum values by so there aren't zeros
#'
#' @import RTMB
#' @returns Returns likelihood values from a logistic normal
#' @keywords internal
#'
dlogistnormal = function(obs, pred, Sigma, perc = 0.85, give_log = TRUE) {

  # Dealing with zeros
  if(any(obs == 0)) {
    # normalize by adding 85% of the minimum observed value in a given year
    obs = (obs + min(obs[obs != 0]) * perc) / sum(obs + min(obs[obs != 0]) * perc)
    pred = (pred + min(obs[obs != 0]) * perc) / sum(pred + min(obs[obs != 0]) * perc)
    # do logistic transformation on observed values
    tmp_Obs = log(obs[-length(obs)])
    tmp_Obs = tmp_Obs - log(obs[length(obs)])
    # do logistic transformation on expected values
    mu = log(pred[-length(pred)]) # remove last bin since it's known
    mu = mu - log(pred[length(pred)]) # calculate log ratio
  } else {
    obs = obs / sum(obs)
    # do logistic transformation on observed values
    tmp_Obs = log(obs[-length(obs)])
    tmp_Obs = tmp_Obs - log(obs[length(obs)])
    # do logistic transformation on expected values
    mu = log(pred[-length(pred)]) # remove last bin since it's known
    mu = mu - log(pred[length(pred)]) # calculate log ratio
  } # if we don't have any zeros
  res = RTMB::dmvnorm(as.vector(tmp_Obs), as.vector(mu), Sigma = Sigma, log = give_log) # fit multivariate normal on log ratios
  return(res)
}


#' Negative binomial that can take non-integer values
#'
#' @param x observations
#' @param give_log whether to give log
#' @param log_mu log mu
#' @param log_var_minus_mu log var minus mu - reparameterize negbin
#'
#' @returns Returns likelihood values from a robust negative binomial
#' @keywords internal
#'
dnbinom_robust_noint <- function(x, log_mu, log_var_minus_mu, give_log = TRUE) {
  mu = exp(log_mu)
  var_minus_mu = exp(log_var_minus_mu)
  k = mu^2 / var_minus_mu # get overdispersion
  logres = lgamma(k+x)-lgamma(k)-lgamma(x+1)+k*log(k)-k*log(mu+k)+x*log(mu)-x*log(mu+k)
  if(give_log) return(logres) else return(exp(logres))
}


#' Poisson that can take non-integer values
#'
#' @param x observations
#' @param pred predicted
#' @param give_log if giving log likelihood
#'
#' @returns Returns likelihood values from a poisson
#' @keywords internal
#'
dpois_noint <- function(x, pred, give_log = TRUE) {
  logres <- -pred + x*log(pred) - lgamma(x+1)
  if(give_log == TRUE) return(logres) else return(exp(logres))
}

#' Get scaled alpha beta parameters for a scaled beta distribution (for steepness) from https://stackoverflow.com/questions/75165770/beta-distribution-with-bounds-at-0-1-0-5
#'
#' @param low lower bound
#' @param high upper bound
#' @param sigma sigma in normal space
#' @param mu mean in normal space
#'
#' @returns Returns beta distribution parmeters with bounds
#' @keywords internal
#'
get_beta_scaled_pars <- function(low, high, mu, sigma) {
  # convert mean and sd to alpha and beta
  scale = high - low
  mean = (mu - low) / scale
  var = (sigma / scale) ^ 2
  # stopifnot(var < mean * (1 - mean))
  var = mean * (1 - mean) / var - 1
  # alpha and beta conversion
  a = mean * var
  b = (1 - mean) * var
  return(c(a,b,low,scale))
}
