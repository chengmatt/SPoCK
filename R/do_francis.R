#' Function to compute francis weights, which is used internally by do_francis_reweighting
#'
#' @param n_regions Number of regions
#' @param n_sexes Number of sexes
#' @param n_fleets Number of fleets (fishery or survey)
#' @param Use Array from data list that specifies whether to use data that year
#' @param ISS Input sample size array
#' @param Pred_array Predicted values array dimensioned by n_regions, n_years, n_ages, n_sexes, n_fleets
#' @param Obs_array Observed values array dimensioned by n_regions, n_years, n_ages, n_sexes, n_fleets
#' @param bins Vector of bins used (age or length)
#' @param weights Array of francis weights (NAs) to apply dimensioned by n_regions, n_years, n_sexes, n_fleets
#' @param comp_type Matrix of composition structure types dimensioned by year and fleet
#'
#' @returns Value for calculated francis weight
#' @export get_francis_weights
#'
#' @examples
#' \dontrun{
#' # Function is used within do_francis_reweighting
#' get_francis_weights(n_regions, n_sexes, n_fleets, Use, ISS, Pred_arry, Obs_array, weights, bins, comp_type)
#' }
get_francis_weights <- function(n_regions,
                                n_sexes,
                                n_fleets,
                                Use,
                                ISS,
                                Pred_array,
                                Obs_array,
                                weights,
                                bins,
                                comp_type
) {

  for(f in 1:n_fleets) {

    data_yrs <- which(Use[,,f] == 1) # get years with data for a given fleet

    # Set up reweighting vectors
    exp_bar <- array(NA, dim = c(n_regions, length(data_yrs), n_sexes)) # mean expected
    obs_bar <- array(NA, dim = c(n_regions, length(data_yrs), n_sexes))  # mean observed
    v_y <- array(NA, dim = c(n_regions, length(data_yrs), n_sexes))  # variance
    w_denom <- array(NA, dim = c(n_regions, length(data_yrs), n_sexes))  # weight factor in denominator

    for(y in data_yrs) {

      # Extract out temporary variables
      tmp_iss_obs <- ISS[,y,,f, drop = FALSE] # temporary ISS
      tmp_exp <- Pred_array[,y,,,f, drop = FALSE] # temporary expected values
      tmp_obs <- Obs_array[,y,,,f, drop = FALSE] # temporary observed values
      yr_alt_idx <- which(data_yrs == y) # get indexing to start from 1

      if(comp_type[y,f] == 3) { # if joint, reorder array
        tmp_exp_jntall <- aperm(tmp_exp, perm = c(3,4,1,2,5)) # reorder array, by ages, sexes, region, year, and fleet
        tmp_obs_jntall <- aperm(tmp_obs, perm = c(3,4,1,2,5)) # reorder array, by ages, sexes, region, year, and fleet
      }

      # If compositions are aggregated across regions and sexes
      if(comp_type[y,f] == 0) {
        exp_bar[1,yr_alt_idx,1] <- sum(bins * as.vector(tmp_exp[1,1,,1,1])) # get mean pred comps
        obs_bar[1,yr_alt_idx,1] <- sum(bins * as.vector(tmp_obs[1,1,,1,1])) # get mean obs comps
        v_y[1,yr_alt_idx,1] <- sum(bins^2*tmp_exp[1,1,,1,1])-exp_bar[1,yr_alt_idx,1]^2 # get variance
        w_denom[1,yr_alt_idx,1] <- (obs_bar[1,yr_alt_idx,1]-exp_bar[1,yr_alt_idx,1])/sqrt(v_y[1,yr_alt_idx,1]/tmp_iss_obs[1,1,1,1]) # get weights
      } # end if aggregated

      # If compositions are split by region and sex
      if(comp_type[y,f] == 1) {
        for(r in 1:n_regions) {
          for(s in 1:n_sexes) {
            exp_bar[r,yr_alt_idx,s] <- sum(bins * as.vector(tmp_exp[r,1,,s,1])) # get mean pred comps
            obs_bar[r,yr_alt_idx,s] <- sum(bins * as.vector(tmp_obs[r,1,,s,1])) # get mean obs comps
            v_y[r,yr_alt_idx,s] <- sum(bins^2*tmp_exp[r,1,,s,1])-exp_bar[r,yr_alt_idx,s]^2 # get variance
            w_denom[r,yr_alt_idx,s] <- (obs_bar[r,yr_alt_idx,s]-exp_bar[r,yr_alt_idx,s])/sqrt(v_y[r,yr_alt_idx,s]/tmp_iss_obs[r,1,s,1]) # get weights
          } # end s loop
        } # end r loop
      } # end if split by region and sex

      # If compositions are split by region, but joint by sex
      if(comp_type[y,f] == 2) {
        for(r in 1:n_regions) {
          exp_bar[r,yr_alt_idx,1] <- sum(rep(bins, n_sexes) * as.vector(tmp_exp[r,1,,,1])) # get mean pred comps
          obs_bar[r,yr_alt_idx,1] <- sum(rep(bins, n_sexes) * as.vector(tmp_obs[r,1,,,1])) # get mean obs comps
          v_y[r,yr_alt_idx,1] <- sum(rep(bins, n_sexes)^2*as.vector(tmp_exp[r,1,,,1]))-exp_bar[r,yr_alt_idx,1]^2 # get variance
          w_denom[r,yr_alt_idx,1] <- (obs_bar[r,yr_alt_idx,1]-exp_bar[r,yr_alt_idx,1])/sqrt(v_y[r,yr_alt_idx,1]/tmp_iss_obs[r,1,1,1]) # get weights
        } # end r loop
      } # end if split by region, joint by sex

      # If compostions are joint by region and sex
      if(comp_type[y,f] == 3) {
        exp_bar[1,yr_alt_idx,1] <- sum(rep(bins, n_regions * n_sexes) * as.vector(tmp_exp_jntall[,,,1,1])) # get mean pred comps
        obs_bar[1,yr_alt_idx,1] <- sum(rep(bins, n_regions * n_sexes) * as.vector(tmp_obs_jntall[,,,1,1])) # get mean obs comps
        v_y[1,yr_alt_idx,1] <- sum(rep(bins, n_regions * n_sexes)^2*as.vector(tmp_exp_jntall[,1,,,1]))-exp_bar[1,yr_alt_idx,1]^2 # get variance
        w_denom[1,yr_alt_idx,1] <- (obs_bar[1,yr_alt_idx,1]-exp_bar[1,yr_alt_idx,1])/sqrt(v_y[1,yr_alt_idx,1]/tmp_iss_obs[1,1,1,1]) # get weights
      } # end if joint by region and sex
    } # end y loop

    # get unique composition types
    unique_comp_type <- unique(comp_type[,f])

    for(j in 1:length(unique_comp_type)) {

      # get year pointer index to subset w_denom and calculate weights separately for each composition type
      year_pointer <- which(comp_type[data_yrs,f] == unique_comp_type[j])

      # if aggregated or joint by region and sex
      if(unique_comp_type[j] %in% c(0,3)) weights[1,data_yrs,1,f] <- 1 / var(w_denom[1,year_pointer,1])

      # if split by sex and region
      if(unique_comp_type[j] == 1) for(r in 1:n_regions) for(s in 1:n_sexes) weights[r,data_yrs,s,f] <- 1 / var(w_denom[r,year_pointer,s])

      # if split by region, joint by sex
      if(unique_comp_type[j] == 2) for(r in 1:n_regions) weights[r,data_yrs,1,f] <- 1 / var(w_denom[r,year_pointer,1])

    } # end j loop
  } # end f loop

  return(weights)
} # end function

#' Function to run Francis reweighting
#'
#' @param rep Report file list
#' @param age_labels Age labels
#' @param len_labels Length labels
#' @param year_labels Year labels
#' @param data List of data inputs
#'
#' @return A list object of francis weights (note that it will be NAs for some, if using jnt composition approaches - i.e., only uses one dimension)
#' @export do_francis_reweighting
#'
#' @examples
#' \dontrun{
#' for(j in 1:5) {
#'
#'   if(j == 1) { # reset weights at 1
#'     data$Wt_FishAgeComps[] <- 1
#'     data$Wt_FishLenComps[] <- 1
#'     data$Wt_SrvAgeComps[] <- 1
#'     data$Wt_SrvLenComps[] <- 1
#'   } else {
#'     data$Wt_FishAgeComps[] <- wts$new_fish_age_wts
#'     data$Wt_FishLenComps[] <- wts$new_fish_len_wts
#'     data$Wt_SrvAgeComps[] <- wts$new_srv_age_wts
#'     data$Wt_SrvLenComps[] <- wts$new_srv_len_wts
#'   }
#'
#'   # make AD model function
#'   sabie_rtmb_model <- RTMB::MakeADFun(cmb(SPoCK_rtmb, data), parameters = parameters, map = mapping, silent = T)
#'
#'   # Now, optimize the function
#'   sabie_optim <- stats::nlminb(sabie_rtmb_model$par, sabie_rtmb_model$fn, sabie_rtmb_model$gr,
#'                                control = list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15))
#'   # newton steps
#'   try_improve <- tryCatch(expr =
#'                             for(i in 1:3) {
#'                               g = as.numeric(sabie_rtmb_model$gr(sabie_optim$par))
#'                               h = optimHess(sabie_optim$par, fn = sabie_rtmb_model$fn, gr = sabie_rtmb_model$gr)
#'                               sabie_optim$par = sabie_optim$par - solve(h,g)
#'                               sabie_optim$objective = sabie_rtmb_model$fn(sabie_optim$par)
#'                             }
#'                           , error = function(e){e}, warning = function(w){w})
#'
#'   rep <- sabie_rtmb_model$report(sabie_rtmb_model$env$last.par.best) # Get report
#'   wts <- do_francis_reweighting(data = data, rep = rep, age_labels = 2:31, len_labels = seq(41, 99, 2), year_labels = 1960:2024)
#' }
#' }
do_francis_reweighting <- function(data,
                                   rep,
                                   age_labels,
                                   len_labels,
                                   year_labels
                                   ) {

  # Get indexing
  n_regions <- data$n_regions
  n_fish_fleets <- data$n_fish_fleets
  n_srv_fleets <- data$n_srv_fleets
  n_sexes <- data$n_sexes

  # Get composition proportions
  comp_prop <- get_comp_prop(data = data,
                             rep = rep,
                             age_labels = age_labels,
                             len_labels = len_labels,
                             year_labels = year_labels)

  # Fishery Ages ------------------------------------------------------------

  # Extract variables
  tmp_Use <- data$UseFishAgeComps
  tmp_ISS <- data$ISS_FishAgeComps
  tmp_Pred <- comp_prop$Pred_FishAge_mat
  tmp_Obs <- comp_prop$Obs_FishAge_mat
  tmp_comp_type <- data$FishAgeComps_Type

  new_fish_age_wts <- data$Wt_FishAgeComps # matrix to store new weights
  new_fish_age_wts[] <- NA

  # Get francis weights here
  new_fish_age_wts[] <- get_francis_weights(
    n_regions = n_regions,
    n_sexes = n_sexes,
    n_fleets = n_fish_fleets,
    Use = tmp_Use,
    ISS = tmp_ISS,
    Pred_array = tmp_Pred,
    Obs_array = tmp_Obs,
    weights = new_fish_age_wts,
    bins = age_labels,
    comp_type = tmp_comp_type
  )

  # Fishery Lengths ------------------------------------------------------------

  # Extract variables
  tmp_Use <- data$UseFishLenComps
  tmp_ISS <- data$ISS_FishLenComps
  tmp_Pred <- comp_prop$Pred_FishLen_mat
  tmp_Obs <- comp_prop$Obs_FishLen_mat
  tmp_comp_type <- data$FishLenComps_Type

  new_fish_len_wts <- data$Wt_FishLenComps # matrix to store new weights
  new_fish_len_wts[] <- NA

  # Get francis weights here
  new_fish_len_wts[] <- get_francis_weights(
    n_regions = n_regions,
    n_sexes = n_sexes,
    n_fleets = n_fish_fleets,
    Use = tmp_Use,
    ISS = tmp_ISS,
    Pred_array = tmp_Pred,
    Obs_array = tmp_Obs,
    weights = new_fish_len_wts,
    bins = len_labels,
    comp_type = tmp_comp_type
  )

  # Survey Ages ------------------------------------------------------------

  # Extract variables
  tmp_Use <- data$UseSrvAgeComps
  tmp_ISS <- data$ISS_SrvAgeComps
  tmp_Pred <- comp_prop$Pred_SrvAge_mat
  tmp_Obs <- comp_prop$Obs_SrvAge_mat
  tmp_comp_type <- data$SrvAgeComps_Type

  new_srv_age_wts <- data$Wt_SrvAgeComps # matrix to store new weights
  new_srv_age_wts[] <- NA

  # Get francis weights here
  new_srv_age_wts[] <- get_francis_weights(
    n_regions = n_regions,
    n_sexes = n_sexes,
    n_fleets = n_srv_fleets,
    Use = tmp_Use,
    ISS = tmp_ISS,
    Pred_array = tmp_Pred,
    Obs_array = tmp_Obs,
    weights = new_srv_age_wts,
    bins = age_labels,
    comp_type = tmp_comp_type
  )

  # Survey Lengths ------------------------------------------------------------

  # Extract variables
  tmp_Use <- data$UseSrvLenComps
  tmp_ISS <- data$ISS_SrvLenComps
  tmp_Pred <- comp_prop$Pred_SrvLen_mat
  tmp_Obs <- comp_prop$Obs_SrvLen_mat
  tmp_comp_type <- data$SrvLenComps_Type

  new_srv_len_wts <- data$Wt_SrvLenComps # matrix to store new weights
  new_srv_len_wts[] <- NA

  # Get francis weights here
  new_srv_len_wts[] <- get_francis_weights(
    n_regions = n_regions,
    n_sexes = n_sexes,
    n_fleets = n_srv_fleets,
    Use = tmp_Use,
    ISS = tmp_ISS,
    Pred_array = tmp_Pred,
    Obs_array = tmp_Obs,
    weights = new_srv_len_wts,
    bins = len_labels,
    comp_type = tmp_comp_type
  )

  return(list(new_fish_age_wts = new_fish_age_wts, new_fish_len_wts = new_fish_len_wts,
              new_srv_age_wts = new_srv_age_wts, new_srv_len_wts = new_srv_len_wts))

} # end function
