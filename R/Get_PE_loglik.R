#' Get selectivity process error log likelihoods (positive)
#'
#' @param PE_model Process error model values (1, 2, 3, 4, and 5) (iid, random walk, 3d marginal, 3d conditional, and 2dar1)
#' @param PE_pars Process error parameters
#' @param ln_devs Deviations
#' @param map_sel_devs selectivity deviations to share
#' @param sel_vals Selectivity values
#'
#' @returns numeric value of log likelihood (in positive space)
#' @keywords internal
#' @import RTMB
Get_sel_PE_loglik <- function(PE_model,
                              PE_pars,
                              ln_devs,
                              map_sel_devs,
                              sel_vals
                              ) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  # Note that the likelihood calculations are positive within the function,
  # because it gets converted to negative outside the wrapper function

  ll = 0 # initialize likelihood

  # find unique selectivity deviations to penalize (sort drops NAs)
  unique_sel_devs = sort(unique(as.vector(map_sel_devs)))

  n_yrs = dim(map_sel_devs)[2] # get years for indexing
  n_ages = dim(map_sel_devs)[3] # get ages for indexing
  n_sexes = dim(map_sel_devs)[4] # get sexes for indexing

  if(PE_model %in% c(1, 2)) {

    for(dev_idx in 1:length(unique_sel_devs)) {

      # figure out where unique sel devs first occur
      idx = which(map_sel_devs == unique_sel_devs[dev_idx], arr.ind = TRUE)[1,]
      y = idx[2] # get unique year index
      i = idx[3] # get unique age or parmeter index
      s = idx[4] # get unique sex index

      # get minimum year process error should be applied to
      min_yr = min(which(map_sel_devs[1,,i,s,1] %in% unique_sel_devs))

      if(PE_model == 1) {
        if(y >= min_yr) ll = ll + RTMB::dnorm(ln_devs[1,y,i,s,1], 0, exp(PE_pars[1,i,s,1]), TRUE)
      } # iid process error

      if(PE_model == 2) {
        if(y == min_yr) ll = ll + RTMB::dnorm(ln_devs[1,y,i,s,1], 0, 5, TRUE) # initialize w/ big value
        else if(y > min_yr) ll = ll + RTMB::dnorm(ln_devs[1,y,i,s,1], ln_devs[1,y-1,i,s,1], exp(PE_pars[1,i,s,1]), TRUE)
      } # end random walk process error

    } # end dev_idx loop

    # Temporal Stability Penalty
    for(y in (min_yr + 1):n_yrs) {
      for(s in 1:n_sexes) {
        for(a in 1:n_ages) {
          year_penalty = log(sel_vals[1,y,a,s,1]) - log(sel_vals[1,y-1,a,s,1])
          ll = ll - year_penalty^2
        } # end a loop
      } # end s loop
    } # end y loop

  } # end iid or random walk process error

  if(PE_model %in% c(3,4,5)) {

    if(PE_model == 3) Var_Type = 0 # marginal variance
    if(PE_model == 4) Var_Type = 1 # conditional variance

    # get first unique combination
    unique_comb = which(map_sel_devs == unique_sel_devs[1], arr.ind = TRUE)[1,]
    # cbind to get all unique combinations here (cbinding first one, so loop starts at 2)
    for(i in 2:length(unique_sel_devs)) unique_comb = cbind(unique_comb, which(map_sel_devs == unique_sel_devs[i], arr.ind = TRUE)[1,])

    # Next, get unique sex deviations
    unique_s = unique(unique_comb[4,])

    for(idx in 1:length(unique_s)) {

      s = unique_s[idx] # get sex index

      # get minimum year process error should be applied to
      min_yr = min(which(map_sel_devs[1,,,s,1] %in% unique_sel_devs))

      # Construct precision matrix for 3d gmrf
      if(PE_model %in% c(3,4)) {
        Q = Get_3d_precision(n_ages = n_ages, # number of ages
                             n_yrs = n_yrs,  # number of years
                             pcorr_age = PE_pars[1,1,s,1], # unconstrained partial correaltion by age
                             pcorr_year = PE_pars[1,2,s,1], # unconstrained partial correaltion by year
                             pcorr_cohort = PE_pars[1,3,s,1], # unconstrained partial correaltion by cohort
                             ln_var_value = PE_pars[1,4,s,1], # log variance
                             Var_Type = Var_Type) # variance type, == 0 (marginal), == 1 (conditional)

        # apply gmrf likelihood
        eps_ay = as.vector(t(ln_devs[1,,,s,1])) # convert to vector
        ll = ll + RTMB::dgmrf(x = eps_ay, mu = 0, Q = Q, log = TRUE)
      } # end if

      # 2dar1 model
      if(PE_model == 5) {
        # Function to constrain values between -1 and 1
        rho_trans = function(x) 2/(1+ exp(-2 * x)) - 1

        # Extract out varaibles and transform into appropriate space
        eps_ya = ln_devs[1,,,s,1] # needs to be in matrix format for dseparable
        rho_a = rho_trans(PE_pars[1,1,s,1]) # correlation across ages
        rho_y = rho_trans(PE_pars[1,2,s,1]) # correlation across years
        sigma2 = exp(PE_pars[1,4,s,1])^2 # get sigma

        # Define 2d scale
        scale = sqrt(sigma2) / sqrt(1 - rho_y^2) / sqrt(1 - rho_a^2)

        # Define ar1 separable functions
        f1 = function(x) RTMB::dautoreg(x, mu = 0, phi = rho_y, log = TRUE)
        f2 = function(x) RTMB::dautoreg(x, mu = 0, phi = rho_a, log = TRUE)
        ll = ll + RTMB::dseparable(f1, f2)(eps_ya, scale = scale)
      } # end if
    } # end idx loop

    # Regularity on ages
    for(s in 1:n_sexes) {
      for (y in min_yr:n_yrs) {
        if (n_ages >= 3) {
          for (a in 2:(n_ages - 1)) {
            age_curvature = log(sel_vals[1,y,a+1,s,1]) - 2 * log(sel_vals[1,y,a,s,1]) + log(sel_vals[1,y,a-1,s,1])
            ll = ll - age_curvature^2
          } # end a loop
        } # end if
      } # end y loop
    } # end s loop

    # Regularity on years
    for(s in 1:n_sexes) {
      for(a in 1:n_ages) {
        if(n_yrs >= 3) {
          for(y in (min_yr+1):(n_yrs - 1)) {
            year_curvature = log(sel_vals[1,y+1,a,s,1]) - 2 * log(sel_vals[1,y,a,s,1]) + log(sel_vals[1,y-1,a,s,1])
            ll = ll - year_curvature^2
          } # end y loop
        } # end if n_yrs >= 3
      } # end a loop
    }

  } # end 3dgrmf or 2dar1 process error

  return(ll)
} # return log likelihood

#' Title Get Movement Process Error Likelihoods
#'
#' @param PE_model Process error model values
#' @param PE_pars Process error parameters
#' @param logit_devs Deviations
#' @param map_move_devs movement deviations to share
#' @param do_recruits_move Whether recruits move (0, don't move, 1 move)
#'
#' @returns numeric value of log likelihood (in positive space)
#' @keywords internal
#' @import RTMB
Get_move_PE_loglik <- function(PE_model,
                               PE_pars,
                               logit_devs,
                               map_move_devs,
                               do_recruits_move
                               ) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  # Note that the likelihood calculations are positive within the function,
  # because it gets converted to negative outside the wrapper function

  ll = 0 # initialize likelihood

  # Get dimensions for curvature penalty
  n_regions_from = dim(map_move_devs)[1]
  n_regions_to = dim(map_move_devs)[2]
  n_yrs = dim(map_move_devs)[3]
  n_ages = dim(map_move_devs)[4]
  n_sexes = dim(map_move_devs)[5]

  # whether recruits move
  age_start = ifelse(do_recruits_move == 0 && n_ages >= 2, 2, 1)

  # Penalize Deviations
  for(rr in 1:n_regions_to) {
    for(r in 1:n_regions_from) {

      if(PE_model == 1) {
        for(y in 1:n_yrs) {
          ll = ll + RTMB::dnorm(logit_devs[r,rr,y,1,1], 0, exp(PE_pars[r,1,1]), TRUE)
        } # end y loop
      } # iid_y

      if(PE_model == 2) {
        for(a in age_start:n_ages) {
          ll = ll + RTMB::dnorm(logit_devs[r,rr,1,a,1], 0, exp(PE_pars[r,a,1]), TRUE)
        } # end a loop
      } # iid_a

      if(PE_model == 3) {
        for(y in 1:n_yrs) {
          for(a in age_start:n_ages) {
            ll = ll + RTMB::dnorm(logit_devs[r,rr,y,a,1], 0, exp(PE_pars[r,a,1]), TRUE)
          } # end a loop
        } # end y loop
      } # iid_y_a

      if(PE_model == 4) {
        for(y in 1:n_yrs) {
          for(s in 1:n_sexes) {
            ll = ll + RTMB::dnorm(logit_devs[r,rr,y,1,s], 0, exp(PE_pars[r,1,s]), TRUE)
          } # end s loop
        } # end y loop
      } # iid_y_s

      if(PE_model == 5) {
        for(a in age_start:n_ages) {
          for(s in 1:n_sexes) {
            ll = ll + RTMB::dnorm(logit_devs[r,rr,1,a,s], 0, exp(PE_pars[r,a,s]), TRUE)
          } # end s loop
        } # end a loop
      } # iid_a_s

      if(PE_model == 6) {
        for(y in 1:n_yrs) {
          for(a in age_start:n_ages) {
            for(s in 1:n_sexes) {
              ll = ll + RTMB::dnorm(logit_devs[r,rr,y,a,s], 0, exp(PE_pars[r,a,s]), TRUE)
            } # end s loop
          } # end a loop
        } # end y loop
      } # iid_y_a_s

    } # end r loop
  } # end rr loop

  return(ll)
}
