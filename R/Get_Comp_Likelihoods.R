#' Gives negative log liklelihood values for composition data for a given year and a given fleet (fishery or survey)
#'
#' @param Exp Expected values (catch at age or survey index at age) indexed for a given year and fleet (structured as a matrix by age and sex)
#' @param Obs Observed values (catch at age or survey index at age) indexed for a given year and fleet (structured as a matrix by age and sex)
#' @param ISS Input sample size indexed for a given year and fleet (structured as a vector w/ sexes)
#' @param Wt_Mltnml Mutlinomial weight (if any) for a given fleet (structured as a vector w/ sexes)
#' @param Comp_Type Composition Parameterization Type (== 0, aggregated comps by sex, == 1, split comps by sex and region (no implicit sex and region ratio information),
#' == 2, joint comps across sexes but split by region (implicit sex ratio information, but not region information))
#' @param Likelihood_Type Composition Likelihood Type (== 0, Multinomial, == 1 Dirichlet Multinomial)
#' @param n_sexes Number of sexes modeled
#' @param age_or_len Age or length comps (== 0, Age, == 1, Length)
#' @param AgeingError Ageing Error matrix
#' @param ln_theta Log theta overdispersion for Dirichlet mutlinomial (scalar or vector depending on if 'Split' or 'Joint')
#' @param n_regions number of regions modeled
#' @param use Vector of 0s and 1s corresponding to regions (==0, don't have obs and dont' use, ==1, have obs and use)
#' @param ln_theta_agg Log overdispersion parameter if comp_type == 0, but we want to estsimate either a dirichlet or multinomial
#' @param comp_agg_type How to aggregate data (if aggregating)
#' @param LN_corr_pars Logistic normal correlation parameters (dimensioned by n_regions, n_sexes, and 3 parameters)
#' @param LN_corr_pars_agg Logistic normal correlation parameters if comps are aggregated (just dimensioned by length of 1 value)
#' @param n_model_bins Number of bins used in the model
#' @param n_obs_bins Number of observed composition bins
#'
#' @return Returns negative log likelihood for composition data (age and/or length)
#' @keywords internal
#'
Get_Comp_Likelihoods = function(Exp,
                                Obs,
                                ISS,
                                Wt_Mltnml,
                                ln_theta_agg,
                                ln_theta,
                                LN_corr_pars = 0,
                                LN_corr_pars_agg = 0,
                                Comp_Type,
                                Likelihood_Type,
                                n_regions,
                                n_model_bins,
                                n_obs_bins,
                                n_sexes,
                                age_or_len,
                                AgeingError,
                                use,
                                comp_agg_type
) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  rho_trans = function(x) 2/(1+ exp(-2 * x)) - 1 # constraint between -1 and 1
  comp_nLL = array(0, dim = c(n_regions, n_sexes)) # initialize nLL here
  const = 0.001 # small constant
  # Filter expectation and observations to regions that have observations
  n_regions_obs_use = sum(use == 1) # get number of regions that have observations

  # Making sure things are correctly formatted (and regions are not dropped)
  Obs = array(Obs, dim = c(n_regions, n_obs_bins, n_sexes))
  Exp = array(Exp, dim = c(n_regions, n_obs_bins, n_sexes)) # using n_obs_bins, because non-square ageing error matrix will collapse to n_obs_bins
  ISS = array(ISS, dim = c(n_regions, n_sexes))
  Wt_Mltnml = array(Wt_Mltnml, dim = c(n_regions, n_sexes))
  ln_theta = array(ln_theta, dim = c(n_regions, n_sexes))
  LN_corr_pars = array(LN_corr_pars, dim = c(n_regions, n_sexes, 3))

  # filter regions that have obs
  Obs = Obs[which(use == 1),,,drop = FALSE]
  Exp = Exp[which(use == 1),,,drop = FALSE]

  # Aggregated comps by sex and region
  if(Comp_Type == 0) {
    if(comp_agg_type == 0) { # aggregated age comps are normalized, aggregated, ageing error, and then normalized again
      # Expected Values
      tmp_Exp = Exp / array(data = rep(colSums(matrix(Exp, nrow = n_model_bins)), each = n_model_bins), dim = dim(Exp)) # normalize by sex and region
      tmp_Exp = matrix(rowSums(matrix(tmp_Exp, nrow = n_model_bins)) / (n_sexes * n_regions), nrow = 1) # take average proportions and transpose
    }

    if(comp_agg_type == 1) tmp_Exp = matrix(rowSums(matrix(Exp, nrow = n_model_bins)) / (n_sexes * n_regions), nrow = 1) # age comps are aggregated, ageing error, and the normalized

    # Expected age bins get collapsed to observed age bins if ageing error is non-square
    if(age_or_len == 0) {
      tmp_Exp = tmp_Exp %*% AgeingError # apply ageing error
      tmp_Exp = as.vector((tmp_Exp) / sum(tmp_Exp)) # renormalize
    }

    if(age_or_len == 1) tmp_Exp = as.vector((tmp_Exp) / sum(tmp_Exp)) # renormalize (lengths)

    # Multinomial likelihood
    if(Likelihood_Type == 0) { # Note that this indexes 1 because it's only a single sex and single region
      tmp_Obs = (Obs[1,,1]) / sum(Obs[1,,1]) # Normalize observed values
      ESS = ISS[1,1] * Wt_Mltnml[1,1] # Effective sample size
      comp_nLL[1,1] = -1 * ESS * sum(((tmp_Obs + const) * log(tmp_Exp + const))) # ADMB multinomial likelihood
      comp_nLL[1,1] = comp_nLL[1,1] - -1 * ESS * sum(((tmp_Obs + const) * log(tmp_Obs + const))) # Multinomial offset (subtract offset from actual likelihood)
    } # end if multinomial likelihood

    if(Likelihood_Type == 1) {
      tmp_Obs = (Obs[1,,1]) / sum(Obs[1,,1]) # Normalize observed values
      comp_nLL[1,1] = -1 * ddirmult(obs = tmp_Obs, pred = tmp_Exp, Ntotal = ISS[1,1], ln_theta = ln_theta_agg, TRUE) # Dirichlet Multinomial likelihood
    } # end if dirichlet multinomial

    if(Likelihood_Type == 2) {
      tmp_Obs = (Obs[1,,1]) / sum(Obs[1,,1]) # Normalize observed values
      Sigma = diag(length(tmp_Obs)-1) * exp(ln_theta_agg)^2
      comp_nLL[1,1] = -1 * dlogistnormal(obs = tmp_Obs, pred = tmp_Exp, Sigma = Sigma, TRUE) # Logistic Normal likelihood (iid)
    } # end if logistic normal (iid)

    if(Likelihood_Type == 3) {
      tmp_Obs = (Obs[1,,1]) / sum(Obs[1,,1]) # Normalize observed values
      LN_corr_b = rho_trans(LN_corr_pars_agg) # correlation by age / length
      Sigma = get_AR1_CorrMat(n_obs_bins, LN_corr_b)
      Sigma = Sigma[-nrow(Sigma), -ncol(Sigma)] * exp(ln_theta_agg)^2 # remove last row and column
      comp_nLL[1,1] = -1 * dlogistnormal(obs = tmp_Obs, pred = tmp_Exp, Sigma = Sigma, TRUE) # Logistic Normal likelihood (1dar1)
    } # end if logistic normal (1dar1)

  } # end if aggregated comps across sexes and regions

  # 'Split' comps by sex and region (no implicit sex ratio information)
  if(Comp_Type == 1) {
    for(s in 1:n_sexes) {
      for(r in 1:n_regions_obs_use) {
        # Expected Values
        if(age_or_len == 0) tmp_Exp = ((Exp[r,,s]) / sum(Exp[r,,s])) %*% AgeingError # Normalize temporary variable (ages)
        if(age_or_len == 1 && comp_agg_type == 0) tmp_Exp = (Exp[r,,s]) / sum(Exp[r,,s]) # Length comps are not normalized prior to age length transition
        if(age_or_len == 1 && comp_agg_type == 1) tmp_Exp = (Exp[r,,s]) # Length comps are normalized prior to age length transition

        # Multinomial likelihood
        if(Likelihood_Type == 0) {
          tmp_Obs = (Obs[r,,s]) / sum(Obs[r,,s]) # Normalize observed temporary variable
          ESS = ISS[r,s] * Wt_Mltnml[r,s] # Effective sample size
          comp_nLL[r,s] = -1 * ESS * sum(((tmp_Obs + const) * log(tmp_Exp + const))) # ADMB multinomial likelihood
          comp_nLL[r,s] = comp_nLL[r,s] - -1 * ESS * sum(((tmp_Obs + const) * log(tmp_Obs + const))) # Multinomial offset (subtract offset from actual likelihood)
        } # end if multinomial likelihood

        if(Likelihood_Type == 1) {
          tmp_Obs = (Obs[r,,s] + const) / sum(Obs[r,,s] + const) # Normalize observed temporary variable
          tmp_Exp = (tmp_Exp + const) / sum(tmp_Exp + const) # Normalize expected temporary variable (lgamma can't take 0s)
          comp_nLL[r,s] = -1 * ddirmult(tmp_Obs, tmp_Exp, ISS[r,s], ln_theta[r,s], TRUE) # Dirichlet Multinomial likelihood
        } # end if dirichlet multinomial

        if(Likelihood_Type == 2) {
          tmp_Obs = Obs[r,,s] / sum(Obs[r,,s]) # extract variable and normalize
          Sigma = diag(length(tmp_Obs)-1) * exp(ln_theta[r,s])^2 # construct sigma
          comp_nLL[r,s] = -1 * dlogistnormal(obs = tmp_Obs, pred = tmp_Exp, Sigma = Sigma, TRUE) # Logistic Normal likelihood
        } # end if logistic normal

        if(Likelihood_Type == 3) {
          tmp_Obs = Obs[r,,s] / sum(Obs[r,,s]) # extract variable and normalize
          LN_corr_b = rho_trans(LN_corr_pars[r,s,1]) # correlation by age / length
          Sigma = get_AR1_CorrMat(n_obs_bins, LN_corr_b)
          Sigma = Sigma[-nrow(Sigma), -ncol(Sigma)] * exp(ln_theta[r,s])^2  # remove last row and column
          comp_nLL[r,s] = -1 * dlogistnormal(obs = tmp_Obs, pred = tmp_Exp, Sigma = Sigma, TRUE) # Logistic Normal likelihood (1dar1)
        } # end if logistic normal (1dar1)

      } # end r loop
    } # end s loop
  } # end if 'Split' comps by sex and region

  # Joint by sex, Split by region
  if(Comp_Type == 2) {
    for(r in 1:n_regions_obs_use) {

      # Expected values
      if(age_or_len == 0) { # if ages
        tmp_Exp = t(as.vector((Exp[r,,])/ sum(Exp[r,,]))) %*% kronecker(diag(n_sexes), AgeingError) # apply ageing error
        tmp_Exp = as.vector((tmp_Exp) / sum(tmp_Exp)) # renormalize to make sure sum to 1
      } # if ages
      if(age_or_len == 1) tmp_Exp = as.vector((Exp[r,,]) / sum((Exp[r,,]))) # Normalize temporary variable (lengths)

      # Multinomial likelihood
      if(Likelihood_Type == 0) { # Indexing by r for a given region since it's 'Split' by region and 1 for sex since it's 'Joint' for sex
        tmp_Obs = as.vector((Obs[r,,]) / sum(Obs[r,,])) # Normalize observed temporary variable
        ESS = ISS[r,1] * Wt_Mltnml[r,1] # Effective sample size
        comp_nLL[r,1] = -1 * ESS * sum(((tmp_Obs + const) * log(tmp_Exp + const))) # ADMB multinomial likelihood
        comp_nLL[r,1] = comp_nLL[r,1] - -1 * ESS * sum(((tmp_Obs + const) * log(tmp_Obs + const))) # Multinomial offset (subtract offset from actual likelihood)
      } # end if multinomial likelihood

      if(Likelihood_Type == 1) {
        tmp_Obs = as.vector((Obs[r,,] + const) / sum(Obs[r,,] + const)) # Normalize observed temporary variable
        tmp_Exp = (tmp_Exp + const) / sum(tmp_Exp + const) # normalize temporary expected variable
        comp_nLL[r,1] = -1 * ddirmult(tmp_Obs, tmp_Exp, ISS[r,1], ln_theta[r,1], TRUE) # Dirichlet Multinomial likelihood
      } # end if dirichlet multinomial

      if(Likelihood_Type == 2) {
        tmp_Obs = Obs[r,,] / sum(Obs[r,,]) # extract temporary observed variable and normalize
        Sigma = diag(length(tmp_Obs)-1) * exp(ln_theta[r,1])^2
        comp_nLL[r,1] = -1 * dlogistnormal(obs = tmp_Obs, pred = tmp_Exp, Sigma = Sigma, TRUE) # Logistic Normal likelihood (iid)
      } # end if logistic normal (iid)

      if(Likelihood_Type == 3) {
        tmp_Obs = Obs[r,,] / sum(Obs[r,,]) # extract variable and normalize
        LN_corr_b = rho_trans(LN_corr_pars[r,1,1]) # correlation by age / length
        Sigma = get_AR1_CorrMat(n_obs_bins * n_sexes, LN_corr_b)
        Sigma = Sigma[-nrow(Sigma), -ncol(Sigma)] * exp(ln_theta[r,1])^2  # remove last row and column
        comp_nLL[r,1] = -1 * dlogistnormal(obs = tmp_Obs, pred = tmp_Exp, Sigma = Sigma, TRUE) # Logistic Normal likelihood (1dar1)
      } # end if logistic normal (1dar1)

      if(Likelihood_Type == 4) {
        tmp_Obs = Obs[r,,] / sum(Obs[r,,]) # extract temporary observed variable and normalize
        LN_corr_b = rho_trans(LN_corr_pars[r,1,1]) # correlation by age / length
        LN_corr_s = rho_trans(LN_corr_pars[r,1,2]) # correlation by sex
        Sigma = kronecker(get_AR1_CorrMat(n_obs_bins, LN_corr_b), get_Constant_CorrMat(n_sexes, LN_corr_s))
        Sigma = Sigma[-nrow(Sigma), -ncol(Sigma)] * exp(ln_theta[r,1])^2 # remove last row and column
        comp_nLL[r,1] = -1 * dlogistnormal(obs = tmp_Obs, pred = tmp_Exp, Sigma = Sigma, TRUE) # Logistic Normal likelihood (1dar1 by age, constant corr by sex)
      } # end if logistic normal (1dar1 by age, constant corr by sex)

    } # end r loop
  } # end if 'Joint' comps by sex, but 'Split' by region

  return(comp_nLL) # return negative log likelihood
} # end function
