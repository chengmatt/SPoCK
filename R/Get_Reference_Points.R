#' Get SPR reference points (Single Region)
#'
#' @param pars Parameter List
#' @param data Data List
#'
#' @keywords internal
#' @import RTMB
#'
#' @examples
#' \dontrun{
#' rep <- obj$report(obj$env$last.par.best) # need to have an RTMB object first
#' # Extract out relevant elements
#' n_ages <- length(data$ages) # number of ages
#' n_years <- length(data$years) # number of years
#' data_list$F_fract_flt <- rep$Fmort[1,n_years,] / sum(rep$Fmort[1,n_years,]) # get fleet F fraction to derive population level selectivity
#' data_list$fish_sel <- array(rep$fish_sel[1,n_years,,1,], dim = c(n_ages, data$n_fish_fleets)) # get female selectivity for all fleets
#' data_list$natmort <- rep$natmort[1,n_years,,1] # get female natural mortality
#' data_list$t_spwn <- t_spwn # specified mortality time up until spawning
#' data_list$WAA <- data$WAA[1,n_years,,1] # weight-at-age for females
#' data_list$MatAA <- data$MatAA[1,n_years,,1] # maturity at age for females
#'
#' data_list$SPR_x <- SPR_x # SPR fraction
#'
#' par_list <- list() # set up parameter list
#' par_list$log_F_x <- log(0.1) # F_x starting value
#'
#' # Make adfun object
#' obj <- RTMB::MakeADFun(cmb(single_region_SPR, data_list), parameters = par_list, map = NULL, silent = TRUE)
#' obj$optim <- stats::nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 1e6, eval.max = 1e6, rel.tol = 1e-15))
#' obj$rep <- obj$report(obj$env$last.par.best) # get report
#' }
single_region_SPR <- function(pars,
                              data
                              ) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  RTMB::getAll(pars, data) # get parameters and data

  n_ages = dim(fish_sel)[1] # number of ages

  # exponentitate reference points to "estimate"
  F_x = exp(log_F_x)

  # set up containers
  SB_age = Nspr = array(0, dim = c(2, n_ages)) # 2 slots in rows, for unfished, and fished at F_x

  # Set up the initial recruits
  Nspr[,1] = 1

  # Loop through and decrement recruit
  for(j in 2:(n_ages - 1)) {
    Nspr[1,j] = Nspr[1,j-1] * exp(-1 * natmort[j-1]) # unfished
    Nspr[2,j] = Nspr[2,j-1] * exp(-1 * (natmort[j-1] + sum(F_fract_flt * F_x * fish_sel[j-1,]))) # fished
  }

  # Accumulate plus group
  Nspr[1,n_ages] = Nspr[1,n_ages-1] * exp(-1 * natmort[n_ages-1])/(1-exp(-1*natmort[n_ages])) # unfished
  Nspr[2,n_ages] = Nspr[2,n_ages-1] * exp(-1 * (natmort[n_ages-1] + sum(F_fract_flt * F_x * fish_sel[n_ages-1,])))/
    (1 - exp(-1 * (natmort[n_ages] + sum(F_fract_flt * F_x * fish_sel[n_ages,])))) # fished

  # Convert numbers at age to spawning biomass at age (t_spwn accounts for mortality up until spawning)
  for(j in 1:n_ages) {
    SB_age[1,j] = Nspr[1,j] * WAA[j] * MatAA[j] * exp(-t_spwn * natmort[j]) # unfished
    SB_age[2,j] = Nspr[2,j] * WAA[j] * MatAA[j] * exp(-t_spwn * (natmort[j] + sum(F_fract_flt * F_x * fish_sel[j,]))) # fished
  }

  # Get spawning biomass per recruit to get spawning potential ratio
  SB0 = sum(SB_age[1,])
  SB_F_x = sum(SB_age[2,])
  SPR = SB_F_x / SB0

  # compute objective function to get F_x
  sprpen = 100 * (SPR - SPR_x)^2

  RTMB::REPORT(SB_age)
  RTMB::REPORT(Nspr)
  RTMB::REPORT(SB0)
  RTMB::REPORT(SB_F_x)
  RTMB::REPORT(F_x)

  return(sprpen)
}

#' Title Get Global SPR Reference Points (Spatial)
#'
#' @param pars Parameter List from RTMB
#' @param data Data List from RTMB
#' @keywords internal
#' @import RTMB
#'
#' @examples
#' \dontrun{
#' SPR_x <- 0.4 # spr fraction
#' rep <- obj$report(obj$env$last.par.best) # need to have an RTMB object first
#' # Extract out relevant elements
#' n_ages <- length(data$ages) # number of ages
#' n_years <- length(data$years) # number of years
#' data_list$F_fract_flt <- rep$Fmort[1,n_years,] / sum(rep$Fmort[1,n_years,]) # get fleet F fraction to derive population level selectivity
#' data_list$fish_sel <- array(rep$fish_sel[1,n_years,,1,], dim = c(n_ages, data$n_fish_fleets)) # get female selectivity for all fleets
#' data_list$natmort <- rep$natmort[1,n_years,,1] # get female natural mortality
#' data_list$t_spwn <- t_spwn # specified mortality time up until spawning
#' data_list$WAA <- data$WAA[1,n_years,,1] # weight-at-age for females
#' data_list$MatAA <- data$MatAA[1,n_years,,1] # maturity at age for females
#' data_list$Rec_Prop <- rep$Rec_trans_prop # unfished recruitment by region
#' data_list$Movement <- array(rep$Movement[,,n_years,,1], dim = c(n_regions, n_regions, n_ages)) # Movement
#' data_list$do_recruits_move <- data$do_recruits_move # whether recruits move
#' data_list$SPR_x <- SPR_x # SPR fraction
#' par_list <- list() # set up parameter list
#' par_list$log_F_x <- log(0.1) # F_x starting value
#' # Make adfun object
#' obj <- RTMB::MakeADFun(cmb(global_SPR, data_list), parameters = par_list, map = NULL, silent = TRUE)
#' obj$optim <- stats::nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 1e6, eval.max = 1e6, rel.tol = 1e-15))
#' obj$rep <- obj$report(obj$env$last.par.best) # get report
#' }
global_SPR <- function(pars,
                       data
                       ) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  RTMB::getAll(pars, data) # get parameters and data

  n_regions = dim(fish_sel)[1] # number of regions
  n_model_ages = dim(fish_sel)[2] # number of model ages
  n_ages = n_model_ages * 10 # get number of ages to iterate through for plus group

  # exponentitate reference points to "estimate"
  F_x = exp(log_F_x)

  # set up containers
  SB_age = Nspr = array(0, dim = c(2, n_regions, n_ages)) # 2 slots in rows, for unfished, and fished at F_x

  # Set up the initial recruits
  Nspr[1,,1] = Rec_Prop
  Nspr[2,,1] = Rec_Prop

  # Loop through, apply movement first, then decrement recruit
  for(j in 2:n_ages) {

    # Get age index to use for demographics
    if(j <= n_model_ages) age_idx = j - 1
    else age_idx = n_model_ages

    # Get temporary values
    tmp_unfished = Nspr[1,,j-1]
    tmp_fished = Nspr[2,,j-1]

    # Apply movement
    if(do_recruits_move == 1 || (do_recruits_move == 0 && j > 2)) {
      tmp_unfished = tmp_unfished %*% Movement[,,age_idx]
      tmp_fished = tmp_fished %*% Movement[,,age_idx]
    }

    # decrement recruits after movement and project forward
    Nspr[1,,j] = tmp_unfished * exp(-1 * natmort[,age_idx]) # unfished
    Nspr[2,,j] = tmp_fished * exp(-1 * (natmort[,age_idx] + apply(F_fract_flt * F_x * fish_sel[,age_idx,,drop = FALSE], 1, sum))) # fished

  } # end j loop

  # Convert numbers at age to spawning biomass at age (t_spwn accounts for mortality up until spawning)
  for(j in 1:n_ages) {
    # Get age index to use for demographics
    if(j < n_model_ages) age_idx = j
    else age_idx = n_model_ages
    SB_age[1,,j] = Nspr[1,,j] * WAA[,age_idx] * MatAA[,age_idx] * exp(-t_spwn * natmort[,age_idx]) # unfished
    SB_age[2,,j] = Nspr[2,,j] * WAA[,age_idx] * MatAA[,age_idx] * exp(-t_spwn * (natmort[,age_idx] + apply(F_fract_flt * F_x * fish_sel[,age_idx,,drop = FALSE], 1, sum))) # fished
  } # end j loop

  # Get spawning biomass per recruit to get spawning potential ratio
  SB0 = sum(SB_age[1,,])
  SB_F_x = sum(SB_age[2,,])
  SPR = SB_F_x / SB0

  # compute objective function to get F_x
  sprpen = 100 * (SPR - SPR_x)^2

  RTMB::REPORT(SB_age)
  RTMB::REPORT(Nspr)
  RTMB::REPORT(SB0)
  RTMB::REPORT(SB_F_x)
  RTMB::REPORT(F_x)

  return(sprpen)
}

#' Title Get FMSY from a Beverton-Holt function (Single Region)
#'
#' @param pars Parameter List
#' @param data Data List
#' @keywords internal
#' @import RTMB
#'
#' @examples
#' \dontrun{
#' rep <- obj$report(obj$env$last.par.best) # need to have an RTMB object first
#' data_list <- list() # set up data list
#' # Extract out relevant elements
#' n_ages <- length(data$ages) # number of ages
#' n_years <- length(data$years) # number of years
#' data_list$F_fract_flt <- rep$Fmort[1,n_years,] / sum(rep$Fmort[1,n_years,]) # get fleet F fraction to derive population level selectivity
#' data_list$fish_sel <- array(rep$fish_sel[1,n_years,,1,], dim = c(n_ages, data$n_fish_fleets)) # get female selectivity for all fleets
#' data_list$natmort <- rep$natmort[1,n_years,,1] # get female natural mortality
#' data_list$t_spwn <- 0 # specified mortality time up until spawning
#' data_list$WAA <- data$WAA[1,n_years,,1] # weight-at-age for females
#' data_list$MatAA <- data$MatAA[1,n_years,,1] # maturity at age for females
#' data_list$h <- rep$h_trans # steepness
#' data_list$R0 <- rep$R0 # unfished recruitment
#'
#' par_list <- list() # set up parameter list
#' par_list$log_Fmsy <- log(0.1) # Fmsy starting value
#'
#' # Make adfun object
#' obj <- RTMB::MakeADFun(cmb(single_region_BH_Fmsy, data_list), parameters = par_list, map = NULL, silent = TRUE)
#' obj$optim <- stats::nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 1e6, eval.max = 1e6, rel.tol = 1e-15))
#' obj$rep <- obj$report(obj$env$last.par.best) # get report
#' obj$sdrep <- sdreport(obj)
#' }
single_region_BH_Fmsy <- function(pars,
                                  data) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  RTMB::getAll(pars, data) # get parameters and data

  n_ages = dim(fish_sel)[1] # number of ages

  # exponentitate reference points to "estimate"
  Fmsy = exp(log_Fmsy)

  # set up containers
  SB_age = Nspr = array(0, dim = c(2, n_ages)) # 2 slots in rows, for unfished, and fished at Fmsy
  CAA = rep(0, n_ages) # catch at age

  # Set up the initial recruits
  Nspr[,1] = 1

  # Loop through and decrement recruits
  for(j in 2:(n_ages - 1)) {
    Nspr[1,j] = Nspr[1,j-1] * exp(-1 * natmort[j-1]) # unfished
    Nspr[2,j] = Nspr[2,j-1] * exp(-1 * (natmort[j-1] + sum(F_fract_flt * Fmsy * fish_sel[j-1,]))) # fished
  }

  # Accumulate plus group
  Nspr[1,n_ages] = Nspr[1,n_ages-1] * exp(-1 * natmort[n_ages-1])/(1-exp(-1*natmort[n_ages])) # unfished
  Nspr[2,n_ages] = Nspr[2,n_ages-1] * exp(-1 * (natmort[n_ages-1] + sum(F_fract_flt * Fmsy * fish_sel[n_ages-1,])))/
    (1 - exp(-1 * (natmort[n_ages] + sum(F_fract_flt * Fmsy * fish_sel[n_ages,])))) # fished

  # Derive spawning biomass per recruit and yield per recruit quantities
  for(j in 1:n_ages) {
    tmp_F = sum(F_fract_flt * Fmsy * fish_sel[j,]) # temporary fishing mortality at age
    tmp_Z = tmp_F + natmort[j] # temporary total mortality at age

    # Convert numbers at age to spawning biomass at age (t_spwn accounts for mortality up until spawning)
    SB_age[1,j] = Nspr[1,j] * WAA[j] * MatAA[j] * exp(-t_spwn * natmort[j]) # unfished
    SB_age[2,j] = Nspr[2,j] * WAA[j] * MatAA[j] * exp(-t_spwn * tmp_Z) # fished

    # Get catch at age to derive yield per recruit
    CAA[j] = Nspr[2,j] * (tmp_F / tmp_Z) * (1 - exp(-tmp_Z)) # Baranov's
  }

  # Get spawning biomass per recruit to get spawning biomass per recruit
  SB0 = sum(SB_age[1,])
  SB_F = sum(SB_age[2,])

  # Get equilibrium recruitment
  Req = (4 * h * R0 * SB_F) / (SB0 * (1 - h) + SB_F * (5 * h - 1))

  # Get yield
  Yield = sum(CAA * WAA) * Req

  # Get Bmsy
  Bmsy = SB_F * Req

  # compute objective function to get Fmsy
  obj_fun = -Yield

  RTMB::REPORT(SB_age)
  RTMB::REPORT(Nspr)
  RTMB::REPORT(SB0)
  RTMB::REPORT(SB_F)
  RTMB::REPORT(Fmsy)
  RTMB::REPORT(Yield)
  RTMB::REPORT(Bmsy)
  RTMB::REPORT(Req)

  return(obj_fun)
}

#' Title Get Global FMSY from a Beverton-Holt (Spatial)
#'
#' @param pars Parameter List
#' @param data Data List
#' @keywords internal
#' @import RTMB
#'
#' @examples
#' \dontrun{
#' rep <- obj$report(obj$env$last.par.best) # need to have an RTMB object first
#' data_list <- list() # set up data list
#' # Extract out relevant elements
#' n_ages <- length(data$ages) # number of ages
#' n_years <- length(data$years) # number of years
#' data_list$F_fract_flt <- rep$Fmort[1,n_years,] / sum(rep$Fmort[1,n_years,]) # get fleet F fraction
#' data_list$fish_sel <- array(rep$fish_sel[1,n_years,,1,], dim = c(n_ages, data$n_fish_fleets))
#' data_list$natmort <- rep$natmort[1,n_years,,1]
#' data_list$t_spwn <- 0
#' data_list$WAA <- data$WAA[1,n_years,,1]
#' data_list$MatAA <- data$MatAA[1,n_years,,1]
#' data_list$h <- rep$h_trans
#' data_list$R0 <- rep$R0
#'
#' par_list <- list()
#' par_list$log_Fmsy <- log(0.1)
#'
#' obj <- RTMB::MakeADFun(cmb(single_region_BH_Fmsy, data_list), parameters = par_list, map = NULL, silent = TRUE)
#' obj$optim <- stats::nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 1e6, eval.max = 1e6, rel.tol = 1e-15))
#' obj$rep <- obj$report(obj$env$last.par.best)
#' obj$sdrep <- sdreport(obj)
#' }
global_BH_Fmsy <- function(pars,
                           data) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  RTMB::getAll(pars, data) # get parameters and data

  n_regions = dim(fish_sel)[1] # number of regions
  n_model_ages = dim(fish_sel)[2] # number of model ages
  n_ages = n_model_ages * 10 # get number of ages to iterate through for plus group

  # exponentitate reference points to "estimate"
  Fmsy = exp(log_Fmsy)

  # set up containers
  SB_age = Nspr = array(0, dim = c(2, n_regions, n_ages)) # 2 slots in rows, for unfished, and fished at Fmsy
  CAA = array(0, c(n_regions, n_ages)) # catch at age

  # Set up the initial recruits
  Nspr[1,,1] = Rec_Prop
  Nspr[2,,1] = Rec_Prop

  # Extend WAA by repeating the last age
  WAA_ext <- cbind(WAA, matrix(WAA[, n_model_ages, drop = FALSE], nrow = n_regions, ncol = n_ages - n_model_ages))

  # Loop through, apply movement first, then decrement recruit
  for(j in 2:n_ages) {

    # Get age index to use for demographics
    if(j <= n_model_ages) age_idx = j - 1
    else age_idx = n_model_ages

    # Get temporary values
    tmp_unfished = Nspr[1,,j-1]
    tmp_fished = Nspr[2,,j-1]

    # Apply movement
    if(do_recruits_move == 1 || (do_recruits_move == 0 && j > 2)) {
      tmp_unfished = tmp_unfished %*% Movement[,,age_idx]
      tmp_fished = tmp_fished %*% Movement[,,age_idx]
    }

    # decrement recruits after movement and project forward
    Nspr[1,,j] = tmp_unfished * exp(-1 * natmort[,age_idx]) # unfished
    Nspr[2,,j] = tmp_fished * exp(-1 * (natmort[,age_idx] + apply(F_fract_flt * Fmsy * fish_sel[,age_idx,,drop = FALSE], 1, sum))) # fished

  } # end j loop

  # Derive spawning biomass per recruit and yield per recruit quantities
  for(j in 1:n_ages) {

    # Get age index to use for demographics
    if(j < n_model_ages) age_idx = j
    else age_idx = n_model_ages

    tmp_F = apply(F_fract_flt * Fmsy * fish_sel[,age_idx,,drop = FALSE], 1, sum) # temporary fishing mortality at age
    tmp_Z = tmp_F + natmort[,age_idx] # temporary total mortality at age

    SB_age[1,,j] = Nspr[1,,j] * WAA[,age_idx] * MatAA[,age_idx] * exp(-t_spwn * natmort[,age_idx]) # unfished
    SB_age[2,,j] = Nspr[2,,j] * WAA[,age_idx] * MatAA[,age_idx] * exp(-t_spwn * tmp_Z) # fished

    # Get catch at age to derive yield per recruit
    CAA[,j] = Nspr[2,,j] * (tmp_F / tmp_Z) * (1 - exp(-tmp_Z)) # Baranov's
  } # end j loop

  # Get spawning biomass per recruit to get spawning potential ratio
  SB0 = sum(SB_age[1,,])
  SB_F = sum(SB_age[2,,])
  SPR = SB_F / SB0

  # Get equilibrium recruitment
  Req = (4 * h * R0 * SB_F) / (SB0 * (1 - h) + SB_F * (5 * h - 1))

  # Get yield
  Yield = sum(CAA * WAA_ext) * Req
  Yield_r = rowSums(CAA * WAA_ext) * Req

  # Get Bmsy
  Bmsy = SB_F * Req

  # compute objective function to get Fmsy
  obj_fun = -Yield

  RTMB::REPORT(SB_age)
  RTMB::REPORT(Nspr)
  RTMB::REPORT(SB0)
  RTMB::REPORT(SB_F)
  RTMB::REPORT(Fmsy)
  RTMB::REPORT(Yield)
  RTMB::REPORT(Yield_r)
  RTMB::REPORT(Bmsy)
  RTMB::REPORT(Req)
  RTMB::REPORT(SPR)

  return(obj_fun)
}

#' Wrapper function to get reference points
#'
#' @param data Data list from RTMB
#' @param rep Report list from RTMB
#' @param SPR_x SPR percentage to target
#' @param t_spwn specified mortality time up until spawning
#' @param type Whether this is a "single_region" reference point or "multi_region"
#' @param what What kind of reference point to use:
#' \describe{
#' \item{SPR}{Spawning Potential Ratio for a single region model}
#' \item{independent_SPR}{Spawning Potential Ratio for a multi region model, without movement }
#' \item{global_SPR}{Global Spawning Potential Ratio for a multi region model, with movement }
#' \item{BH_MSY}{MSY reference points derived from a Beverton-Holt for a single region model}
#' \item{independent_BH_MSY}{MSY reference points derived from a Beverton-Holt, without movement }
#' \item{global_BH_MSY}{MSY reference points derived from a Beverton-Holt, with movement }
#' }
#' @param sex_ratio_f Sex ratio for females used to compute biological reference points
#' @param calc_rec_st_yr The first year in which mean recruitment is computed for
#' @param rec_age Actual recruitment age value
#'
#' @importFrom stats nlminb
#' @import RTMB
#'
#' @returns A list object of fishing and biological reference points
#' @export Get_Reference_Points
#'
#' @examples
#' \dontrun{
#' f_40 <- Get_Reference_Points(data = data,
#' rep = rep,
#' SPR_x = 0.4,
#' t_spwn = 0,
#' type = "single_region",
#' what = "SPR")
#' f_40$F_x # estimated reference point
#' }
Get_Reference_Points <- function(data,
                                 rep,
                                 SPR_x = NULL,
                                 t_spwn = 0,
                                 sex_ratio_f = 0.5,
                                 calc_rec_st_yr = 1,
                                 rec_age = 1,
                                 type,
                                 what
                                 ) {

  f_ref_pt <- vector() # set up storage
  b_ref_pt <- vector() # set up storage

  if(type == "single_region") {

    if(!what %in% c("SPR", "BH_MSY")) stop("what is not correctly specified! Should be SPR, BH_MSY for type = single_region")

    data_list <- list() # set up data list
    # Extract out relevant elements
    n_ages <- length(data$ages) # number of ages
    n_years <- length(data$years) # number of years
    data_list$F_fract_flt <- rep$Fmort[1,n_years,] / sum(rep$Fmort[1,n_years,]) # get fleet F fraction to derive population level selectivity
    data_list$fish_sel <- array(rep$fish_sel[1,n_years,,1,], dim = c(n_ages, data$n_fish_fleets)) # get female selectivity for all fleets
    data_list$natmort <- rep$natmort[1,n_years,,1] # get female natural mortality
    data_list$t_spwn <- t_spwn # specified mortality time up until spawning
    data_list$WAA <- data$WAA[1,n_years,,1] # weight-at-age for females
    data_list$MatAA <- data$MatAA[1,n_years,,1] # maturity at age for females

    if(what == 'SPR') {
      data_list$SPR_x <- SPR_x # SPR fraction

      par_list <- list() # set up parameter list
      par_list$log_F_x <- log(0.1) # F_x starting value

      # Make adfun object
      obj <- RTMB::MakeADFun(cmb(single_region_SPR, data_list), parameters = par_list, map = NULL, silent = TRUE)
      obj$optim <- stats::nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 1e6, eval.max = 1e6, rel.tol = 1e-15))
      obj$rep <- obj$report(obj$env$last.par.best) # get report

      # Output reference points
      f_ref_pt[1] <- obj$rep$F_x
      b_ref_pt[1] <- obj$rep$SB_F_x * sex_ratio_f * mean(rep$Rec[1,calc_rec_st_yr:(n_years - rec_age)])

    } # end SPR reference points

    if(what == 'BH_MSY') {

      # extract out beverton-holt parameters
      data_list$h <- rep$h_trans # steepness
      data_list$R0 <- rep$R0 # unfished recruitment

      par_list <- list() # set up parameter list
      par_list$log_Fmsy <- log(0.1) # Fmsy starting value

      # make adfun ect
      obj <- RTMB::MakeADFun(cmb(single_region_BH_Fmsy, data_list), parameters = par_list, map = NULL, silent = TRUE)
      obj$optim <- stats::nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 1e6, eval.max = 1e6, rel.tol = 1e-15))
      obj$rep <- obj$report(obj$env$last.par.best) # get report

      # Output reference points
      f_ref_pt[1] <- obj$rep$Fmsy
      b_ref_pt[1] <- obj$rep$Bmsy * sex_ratio_f
    }
  }

  if(type == 'multi_region') {

    if(!what %in% c("independent_SPR", "independent_BH_MSY",
                    "global_SPR", "global_BH_MSY"))
      stop("what is not correctly specified! Should be independent_SPR, independent_BH_MSY, global_SPR, global_BH_MSY for type = multi_region")

    data_list <- list() # set up data list

    if(what == "independent_SPR") {
      for(r in 1:data$n_regions) {

        # Extract out relevant elements for a given region
        n_years <- length(data$years) # number of years
        n_ages <- length(data$ages) # number of ages
        data_list$F_fract_flt <- rep$Fmort[r,n_years,] / sum(rep$Fmort[r,n_years,]) # get fleet F fraction to derive population level selectivity
        data_list$fish_sel <- array(rep$fish_sel[r,n_years,,1,], dim = c(n_ages, data$n_fish_fleets)) # get female selectivity for all fleets
        data_list$natmort <- rep$natmort[r,n_years,,1] # get female natural mortality
        data_list$t_spwn <- t_spwn # specified mortality time up until spawning
        data_list$WAA <- data$WAA[r,n_years,,1] # weight-at-age for females
        data_list$MatAA <- data$MatAA[r,n_years,,1] # maturity at age for females

        data_list$SPR_x <- SPR_x # SPR fraction

        par_list <- list() # set up parameter list
        par_list$log_F_x <- log(0.1) # F_x starting value

        # Make adfun object
        tmp_obj <- RTMB::MakeADFun(cmb(single_region_SPR, data_list), parameters = par_list, map = NULL, silent = TRUE)
        tmp_obj$optim <- stats::nlminb(tmp_obj$par, tmp_obj$fn, tmp_obj$gr, control = list(iter.max = 1e6, eval.max = 1e6, rel.tol = 1e-15))
        tmp_obj$rep <- tmp_obj$report(tmp_obj$env$last.par.best) # get report

        # Output reference points
        f_ref_pt[r] <- tmp_obj$rep$F_x
        b_ref_pt[r] <- tmp_obj$rep$SB_F_x * sex_ratio_f * mean(rep$Rec[r,calc_rec_st_yr:(n_years - rec_age)])

      } # end r loop
    } # end independent_SPR

    if(what == "independent_BH_MSY") {
      for(r in 1:data$n_regions) {

        # Extract out relevant elements for a given region
        n_years <- length(data$years) # number of years
        n_ages <- length(data$ages) # number of ages
        data_list$F_fract_flt <- rep$Fmort[r,n_years,] / sum(rep$Fmort[r,n_years,]) # get fleet F fraction to derive population level selectivity
        data_list$fish_sel <- array(rep$fish_sel[r,n_years,,1,], dim = c(n_ages, data$n_fish_fleets)) # get female selectivity for all fleets
        data_list$natmort <- rep$natmort[r,n_years,,1] # get female natural mortality
        data_list$t_spwn <- t_spwn # specified mortality time up until spawning
        data_list$WAA <- data$WAA[r,n_years,,1] # weight-at-age for females
        data_list$MatAA <- data$MatAA[r,n_years,,1] # maturity at age for females
        data_list$h <- rep$h_trans[r] # steepness
        data_list$R0 <- rep$R0 * rep$Rec_trans_prop[r] # unfished recruitment by region

        par_list <- list() # set up parameter list
        par_list$log_Fmsy <- log(0.1) # Fmsy starting value

        # Make adfun object
        tmp_obj <- RTMB::MakeADFun(cmb(single_region_BH_Fmsy, data_list), parameters = par_list, map = NULL, silent = TRUE)
        tmp_obj$optim <- stats::nlminb(tmp_obj$par, tmp_obj$fn, tmp_obj$gr, control = list(iter.max = 1e6, eval.max = 1e6, rel.tol = 1e-15))
        tmp_obj$rep <- tmp_obj$report(tmp_obj$env$last.par.best) # get report

        # Output reference points
        f_ref_pt[r] <- tmp_obj$rep$Fmsy
        b_ref_pt[r] <- tmp_obj$rep$Bmsy * sex_ratio_f

      } # end r loop
    } # end independent_SPR

    if(what == 'global_SPR') {

      # Extract out relevant elements for a given region
      n_ages <- length(data$ages) # number of ages to iterate through
      n_years <- length(data$years) # number of years
      n_regions <- data$n_regions # number of regions
      data_list$F_fract_flt <- rep$Fmort[,n_years,,drop = FALSE] / apply(rep$Fmort[,n_years,,drop = FALSE], 1, sum) # get fleet F fraction to derive population level selectivity
      data_list$fish_sel <- array(rep$fish_sel[,n_years,,1,], dim = c(n_regions, n_ages, data$n_fish_fleets)) # get female selectivity for all fleets
      data_list$natmort <- rep$natmort[,n_years,,1] # get female natural mortality
      data_list$t_spwn <- t_spwn # specified mortality time up until spawning
      data_list$WAA <- data$WAA[,n_years,,1] # weight-at-age for females
      data_list$MatAA <- data$MatAA[,n_years,,1] # maturity at age for females
      data_list$Movement <- array(rep$Movement[,,n_years,,1], dim = c(n_regions, n_regions, n_ages)) # Movement
      data_list$do_recruits_move <- data$do_recruits_move # whether recruits move
      data_list$Rec_Prop <- rep$Rec_trans_prop # recruitment proportions

      data_list$SPR_x <- SPR_x # SPR fraction

      par_list <- list() # set up parameter list
      par_list$log_F_x <- log(0.1) # F_x starting value

      # make adfn object
      obj <- RTMB::MakeADFun(cmb(global_SPR, data_list), parameters = par_list, map = NULL, silent = TRUE)
      obj$optim <- stats::nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 1e6, eval.max = 1e6, rel.tol = 1e-15))
      obj$rep <- obj$report(obj$env$last.par.best) # get report

      # output reference points
      f_ref_pt <- rep(obj$rep$F_x, n_regions)
      b_ref_pt <- obj$rep$SB_F_x * sex_ratio_f * rowMeans(rep$Rec[,calc_rec_st_yr:(n_years - rec_age)])

    } # end global SPR

    if(what == 'global_BH_MSY') {

      # extract out dimensions
      n_ages <- length(data$ages) # number of ages to iterate through
      n_years <- length(data$years) # number of years
      n_regions <- data$n_regions # number of regions
      data_list$F_fract_flt <- rep$Fmort[,n_years,,drop = FALSE] / apply(rep$Fmort[,n_years,,drop = FALSE], 1, sum) # get fleet F fraction to derive population level selectivity
      data_list$fish_sel <- array(rep$fish_sel[,n_years,,1,], dim = c(n_regions, n_ages, data$n_fish_fleets)) # get female selectivity for all fleets
      data_list$natmort <- rep$natmort[,n_years,,1] # get female natural mortality
      data_list$t_spwn <- t_spwn # specified mortality time up until spawning
      data_list$WAA <- data$WAA[,n_years,,1] # weight-at-age for females
      data_list$MatAA <- data$MatAA[,n_years,,1] # maturity at age for females
      data_list$Movement <- array(rep$Movement[,,n_years,,1], dim = c(n_regions, n_regions, n_ages)) # Movement
      data_list$do_recruits_move <- data$do_recruits_move # whether recruits move
      data_list$Rec_Prop <- rep$Rec_trans_prop # recruitment proportions
      data_list$h <- mean(rep$h_trans) # steepness
      data_list$R0 <- rep$R0  # unfished recruitment

      par_list <- list() # set up parameter list
      par_list$log_Fmsy <- log(0.1) # Fmsy starting value

      # Make adfun object
      obj <- RTMB::MakeADFun(cmb(global_BH_Fmsy, data_list), parameters = par_list, map = NULL, silent = TRUE)
      obj$optim <- stats::nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 1e6, eval.max = 1e6, rel.tol = 1e-15))
      obj$rep <- obj$report(obj$env$last.par.best) # get report

      # Output reference points
      f_ref_pt <- rep(obj$rep$Fmsy, n_regions)
      b_ref_pt <- obj$rep$Bmsy * sex_ratio_f * rep$Rec_trans_prop
    }

  } # end multi region

  return(list(f_ref_pt = f_ref_pt,
              b_ref_pt = b_ref_pt))
}


