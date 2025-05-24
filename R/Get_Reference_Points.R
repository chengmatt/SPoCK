#' Get SPR reference points
#'
#' @param pars Parameter List
#' @param data Data List
#'
#' @keywords internal
#' @import RTMB
#'
#' @examples
#' \dontrun{
#' # extract out elements from report
#' F_fract_flt <- rep$Fmort[1,65,] / sum(rep$Fmort[1,65,])
#' fish_sel <- rep$fish_sel[1,65,,1,]
#' natmort <- rep$natmort[1,65,,1]
#' t_spwn <- 0
#' WAA <- rtmb_data$WAA[1,65,,1]
#' MatAA <- rtmb_data$MatAA[1,65,,1]
#' SPR_x <- 0.4 # SPR % to target
#' # set up data list
#' data <- list(F_fract_flt = F_fract_flt,
#'              fish_sel = fish_sel,
#'              natmort = natmort,
#'              t_spwn = t_spwn,
#'              WAA = WAA,
#'              MatAA = MatAA,
#'              SPR_x = SPR_x)
#' pars <- list(log_F_x = 0.4)
#' # make adfun and run model
#' spr_ad <- RTMB::MakeADFun(cmb(single_region_SPR, data), parameters = pars, map = NULL)
#' optim <- stats::nlminb(spr_ad$par, spr_ad$fn, spr_ad$gr,
#'                        control = list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15))
#' exp(optim$par) # f_x
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
    SB_age[2,j] = Nspr[2,j] * WAA[j] * MatAA[j] * exp(-t_spwn * (natmort[j] + sum(F_fract_flt * F_x * fish_sel[n_ages-1,]))) # fished
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

#' Wrapper function to get reference points
#'
#' @param data Data list from RTMB
#' @param rep Report list from RTMB
#' @param SPR_x SPR percentage to target
#' @param t_spwn specified mortality time up until spawning
#' @param type Whether this is a "single_region" reference point or "multi_region"
#' @param what What kind of reference point "SPR" and "SPR_independent" (used for multi region, where SPR is calcualted for each region independently)
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
                                 type = "single_region",
                                 what = "SPR"
                                 ) {

  f_ref_pt <- vector() # set up storage
  b_ref_pt <- vector() # set up storage

  if(type == "single_region") {

    data_list <- list() # set up data list
    # Extract out relevant elements
    n_ages <- length(data$ages) # number of ages
    n_years <- length(data$years) # number of years
    data_list$F_fract_flt <- rep$Fmort[1,n_years,] / sum(rep$Fmort[1,n_years,]) # get fleet F fraction to derive population level selectivity
    data_list$fish_sel <- array(rep$fish_sel[1,n_years,,1,], dim = c(n_ages, data$n_fish_fleets)) # get female selectivity for all fleets
    data_list$natmort <- rep$natmort[1,n_years,,1] # get female natural mortality
    data_list$t_spwn <- t_spwn # specified mortality time up until spawning
    data_list$WAA <- data$WAA[1,n_years,,1] # weight-at-age for females
    data_list$MatAA <- data$MatAA[1,n_years,,1] # maturity at age for males

    if(what == 'SPR') {
      data_list$SPR_x <- SPR_x # SPR fraction

      par_list <- list() # set up parameter list
      par_list$log_F_x <- log(0.1) # F_x starting value

      # Make adfun object
      spr_ad <- RTMB::MakeADFun(cmb(single_region_SPR, data_list), parameters = par_list, map = NULL, silent = TRUE)
      spr_ad$optim <- stats::nlminb(spr_ad$par, spr_ad$fn, spr_ad$gr, control = list(iter.max = 1e6, eval.max = 1e6, rel.tol = 1e-15))
      spr_ad$rep <- spr_ad$report(spr_ad$env$last.par.best) # get report

      # Output reference points
      f_ref_pt[1] <- spr_ad$rep$F_x
      b_ref_pt[1] <- spr_ad$rep$SB_F_x * sex_ratio_f * mean(rep$Rec[1,calc_rec_st_yr:(n_years - rec_age)])

    } # end SPR reference points
  }

  if(type == 'multi_region') {

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
        data_list$MatAA <- data$MatAA[r,n_years,,1] # maturity at age for males

        data_list$SPR_x <- SPR_x # SPR fraction

        par_list <- list() # set up parameter list
        par_list$log_F_x <- log(0.1) # F_x starting value

        # Make adfun object
        tmp_spr_ad <- RTMB::MakeADFun(cmb(single_region_SPR, data_list), parameters = par_list, map = NULL, silent = TRUE)
        tmp_spr_ad$optim <- stats::nlminb(tmp_spr_ad$par, tmp_spr_ad$fn, tmp_spr_ad$gr, control = list(iter.max = 1e6, eval.max = 1e6, rel.tol = 1e-15))
        tmp_spr_ad$rep <- tmp_spr_ad$report(tmp_spr_ad$env$last.par.best) # get report

        # Output reference points
        f_ref_pt[r] <- tmp_spr_ad$rep$F_x
        b_ref_pt[r] <- tmp_spr_ad$rep$SB_F_x * sex_ratio_f * mean(rep$Rec[r,calc_rec_st_yr:(n_years - rec_age)])

      } # end r loop
    } # end independent_SPR

  } # end multi region

  return(list(f_ref_pt = f_ref_pt,
              b_ref_pt = b_ref_pt))
}

