#' Get selectivity process error log likelihoods (positive)
#'
#' @param PE_model Process error model values (1, 2, 3, 4, and 5) (iid, random walk, 3d marginal, 3d conditional, and 2dar1)
#' @param PE_pars Process error parameters
#' @param ln_devs Deviations
#' @param map_sel_devs selectivity deviations to share
#'
#' @returns numeric value of log likelihood (in positive space)
#' @keywords internal
#' @import RTMB
Get_sel_PE_loglik <- function(PE_model,
                              PE_pars,
                              ln_devs,
                              map_sel_devs
                              ) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  # Note that the likelihood calculations are positive within the function,
  # because it gets converted to negative outside the wrapper function

  ll = 0 # initialize likelihood

  # find unique selectivity deviations to penalize (sort drops NAs)
  unique_sel_devs = sort(unique(as.vector(map_sel_devs)))

  if(PE_model %in% c(1, 2)) {
    for(dev_idx in 1:length(unique_sel_devs)) {

      # figure out where unique sel devs first occur
      idx = which(map_sel_devs == unique_sel_devs[dev_idx], arr.ind = TRUE)[1,]
      r = idx[1] # get unique region index
      y = idx[2] # get unique year index
      i = idx[3] # get unique age or parmeter index
      s = idx[4] # get unique sex index
      f = idx[5] # get unique fleet index

      if(PE_model == 1) {
        ll = ll + RTMB::dnorm(ln_devs[r,y,i,s,1], 0, exp(PE_pars[r,i,s,1]), TRUE)
      } # iid process error

      if(PE_model == 2) { # random walk
        if(y == 1) {
          ll = ll + RTMB::dnorm(ln_devs[r,1,i,s,1], 0, 50, TRUE) # if y == 1, initialize with large sigma on dnorm
        } else {
          ll = ll + RTMB::dnorm(ln_devs[r,y,i,s,1], ln_devs[r,y-1,i,s,1], exp(PE_pars[r,i,s,1]), TRUE)
        } # end else
      } # end random walk process error

    } # end dev_idx loop
  } # end iid or random walk process error

  if(PE_model %in% c(3,4,5)) {

    if(PE_model == 3) Var_Type = 0 # marginal variance
    if(PE_model == 4) Var_Type = 1 # conditional variance

    # Get indexing for constructor algorithim
    n_yrs = dim(map_sel_devs)[2]
    n_ages = dim(map_sel_devs)[3]

    # get first unique combination
    unique_comb = which(map_sel_devs == unique_sel_devs[1], arr.ind = TRUE)[1,]
    # cbind to get all unique combinations here (cbinding first one, so loop starts at 2)
    for(i in 2:length(unique_sel_devs)) unique_comb = cbind(unique_comb, which(map_sel_devs == unique_sel_devs[i], arr.ind = TRUE)[1,])

    # Next, get unique region, sex combinations
    unique_rs = expand.grid(unique(unique_comb[1,]), unique(unique_comb[4,]))

    for(idx in 1:nrow(unique_rs)) {
      r = unique_rs[idx,1] # get region index
      s = unique_rs[idx,2] # get sex index

      # Construct precision matrix for 3d gmrf
      if(PE_model %in% c(3,4)) {
        Q = Get_3d_precision(n_ages = n_ages, # number of ages
                             n_yrs = n_yrs,  # number of years
                             pcorr_age = PE_pars[r,1,s,1], # unconstrained partial correaltion by age
                             pcorr_year = PE_pars[r,2,s,1], # unconstrained partial correaltion by year
                             pcorr_cohort = PE_pars[r,3,s,1], # unconstrained partial correaltion by cohort
                             ln_var_value = PE_pars[r,4,s,1], # log variance
                             Var_Type = Var_Type) # variance type, == 0 (marginal), == 1 (conditional)

        # apply gmrf likelihood
        eps_ay = as.vector(t(ln_devs[r,,,s,1])) # convert to vector
        ll = ll + RTMB::dgmrf(x = eps_ay, mu = 0, Q = Q, log = TRUE)
      } # end if

      # 2dar1 model
      if(PE_model == 5) {
        # Function to constrain values between -1 and 1
        rho_trans = function(x) 2/(1+ exp(-2 * x)) - 1

        # Extract out varaibles and transform into appropriate space
        eps_ya = ln_devs[r,,,s,1] # needs to be in matrix format for dseparable
        rho_a = rho_trans(PE_pars[r,1,s,1]) # correaltion across ages
        rho_y = rho_trans(PE_pars[r,2,s,1]) # correaltion across ages
        sigma2 = exp(PE_pars[r,4,s,1])^2 # get sigma

        # Define 2d unit variance
        unit_var = sqrt(sigma2) / sqrt(1 - rho_y^2) / sqrt(1 - rho_a^2)

        # Define ar1 separable functions
        f1 = function(x) RTMB::dautoreg(x, mu = 0, phi = rho_y, log = TRUE)
        f2 = function(x) RTMB::dautoreg(x, mu = 0, phi = rho_a, log = TRUE)
        ll = ll + RTMB::dseparable(f1, f2)(eps_ya, scale = unit_var)
      } # end if

    } # end idx loop
  } # end 3dgrmf or 2dar1 process error

  return(ll)
} # return log likelihood
