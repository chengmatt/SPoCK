#' Simulate logistic normal variables
#'
#' @param exp Expected values
#' @param pars Parameters for a logistic normal (iid == 1 parameter, AR1 == 2 parameters, 2D, by age and sex == 3 parameters, 3D, by age, sex, and region == 4 parameters)
#' @param comp_like Likelihood structure (iid == 2, ar1 == 3, 2d == 4, 3d == 5)
#' @importFrom MASS mvrnorm
#' @keywords internal
#'
rlogistnormal <- function(exp,
                          pars,
                          comp_like
                          ) {
  # set up expected value vector
  mu <- log(exp[-length(exp)]) # remove last bin since it's known
  mu <- mu - log(exp[length(exp)]) # calculate log ratio

  # if iid logistic normal
  if(comp_like == 2) {
    Sigma <- diag(length(exp)-1) # set up sigma
    diag(Sigma) <- pars[1]^2 # input parameter
  } # end if iid logistic normal

  # if logistic normal, AR1 by age
  if(comp_like == 3) {
    Sigma <- get_AR1_CorrMat(n_ages, pars[2]) * pars[1]^2
    Sigma <- Sigma[-nrow(Sigma), -ncol(Sigma)] # remove last row and column
  } # end if iid logistic normal

  # if logistic normal, AR1 by age, constant correlation by sex
  if(comp_like == 4) {
    Sigma <- kronecker(get_AR1_CorrMat(n_ages, pars[2]),
                       get_Constant_CorrMat(n_sexes, pars[3])) * pars[1]^2
    Sigma <- Sigma[-nrow(Sigma), -ncol(Sigma)] # remove last row and column
  }

  # if logistic normal, AR1 by age, constant correlation by sex and constant correlation by region
  if(comp_like == 5) {
    Sigma <- kronecker(kronecker(get_AR1_CorrMat(n_ages, pars[2]), get_Constant_CorrMat(n_sexes, pars[3])),
                       get_Constant_CorrMat(n_regions, pars[4])) * pars[1]^2
    Sigma <- Sigma[-nrow(Sigma), -ncol(Sigma)] # remove last row and column
  }

  x <- MASS::mvrnorm(1, mu, Sigma) # simulate from mvnorm (does not sum to 1) and length k
  p <- exp(x)/(1 + sum(exp(x))) # do additive transformation length k and does not sum to 1
  p <- c(p, 1 - sum(p)) # output now so it sums to 1

  return(p)
}


#' Title Simulate dirichlet multinomial draws
#'
#' @param n Number of sims
#' @param N Sum of observations
#' @param alpha Concentration parameter
#' @importFrom stats rgamma rmultinom
#' @keywords internal
rdirM <- function(n, N, alpha) {

  # Get dirichlet draws
  rdirichlet <- function(alpha) {
    x <- stats::rgamma(length(alpha), shape = alpha, scale = 1)
    return(x / sum(x))
  }

  # Generate DM samples
  result <- replicate(n, {
    p <- rdirichlet(alpha)
    counts <- stats::rmultinom(1, size = N, prob = p)
    as.vector(counts)
  })

  return(result)
}

#' Title Generate Recruitment Values based on Inverse Gaussian Distribution
#'
#' @param sims Number of Simulations
#' @param recruitment Recruitment vector
#' @importFrom stats rnorm runif
#' @returns Random variables following an inverse gaussian
#' @keywords internal
rinvgauss_rec <- function(sims,
                          recruitment
                          ) {

  AMeanRec <- mean(recruitment) # get arithmetic mean recruitment
  HMeanRec <- 1 / mean(1 / recruitment) # get harmonic mean recruitment

  # define parameters for inverse gaussian
  gamma <- AMeanRec / HMeanRec
  gi_beta <- AMeanRec
  delta <- 1 / (gamma - 1)
  cvrec <- sqrt(1 / delta)

  # Generate random variables with transformation
  psi <- stats::rnorm(sims,0,1)^2 # generate squared random normal
  omega <- gi_beta * (1 + (psi - sqrt(4 * delta * psi + psi^2)) / (2 * delta))
  zeta <- gi_beta * (1 + (psi + sqrt(4 * delta * psi + psi^2)) / (2 * delta))
  gtheta <- gi_beta / (gi_beta + omega)

  unifs <- stats::runif(sims) # generate uniform
  rv <- ifelse(unifs <= gtheta, omega, zeta) # generate random draws based on mixture

  return(rv)
}

