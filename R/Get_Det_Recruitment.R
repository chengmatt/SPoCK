#' Get Deterministic Recruitment
#'
#' @param recruitment_model == 0, mean recruitment, == 1 beverton holt recruitment with steepness
#' @param R0 virgin or mean recruitment (global)
#' @param h Vector of steepness values by n_regions
#' @param n_ages number of ages
#' @param WAA Weight at age by region and age
#' @param MatAA Maturity by region and age
#' @param natmort Natural moratliaty by region and age
#' @param SSB_vals SSB matrix by regiona and year
#' @param y year
#' @param rec_lag recruitment lag for indexing SSB year
#' @param recruitment_dd Recruitment density dependence (0 == local, 1 == global)
#' @param Rec_Prop Recruitment proportions to allocate global R0 if local density dependence
#' @param n_regions Number of regions
#' @export Get_Det_Recruitment
#' @returns Vector of n_regions of deterministic recruitment values from mean recruitment or beverton holt
#' @keywords internal
#'
Get_Det_Recruitment <- function(recruitment_model,
                                recruitment_dd,
                                y,
                                rec_lag,
                                R0,
                                Rec_Prop,
                                h,
                                n_regions,
                                n_ages,
                                WAA,
                                MatAA,
                                natmort,
                                SSB_vals
                                ) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  if(recruitment_model == 0) rec = R0 * Rec_Prop # mean recruitment apportioned across n_regions

  # Beverton-Holt
  if(recruitment_model == 1) {

    # Storage for recruitment and S0
    rec = S0 = rep(0, n_regions)
    # Calculate unexploited naa per recruit
    tmp_naa = array(0, dim = c(n_regions, n_ages))
    tmp_naa[,1] = 1

    # Restructure some stuff
    tmp_natmort <- array(natmort, dim = c(n_regions, n_ages))
    tmp_WAA <- array(WAA, dim = c(n_regions, n_ages))
    tmp_MatAA <- array(MatAA, dim = c(n_regions, n_ages))

    # Derive virgin spawning biomass
    for(r in 1:n_regions) {
      for(a in 2:n_ages) tmp_naa[r,a] = tmp_naa[r,a-1] * exp(-tmp_natmort[r,a-1])
      if(recruitment_dd == 0) S0[r] = sum(tmp_naa[r,] * tmp_WAA[r,] * tmp_MatAA[r,]) * R0 * Rec_Prop[r] # Get local S0
      if(recruitment_dd == 1) S0[r] = sum(tmp_naa[r,] * tmp_WAA[r,] * tmp_MatAA[r,]) * R0 # Get global S0
    } # end r loop

    # get SSB to use to predict recruitment
    if(y <= rec_lag) SSB = S0 else SSB = SSB_vals[,y-rec_lag]

    # Get recruitment based on SSB and R0
    for(r in 1:n_regions) {

      # Local Density Dependence
      if(recruitment_dd == 0) {
        local_R0 = R0 * Rec_Prop[r] # get local R0 based on recruitment proportions
        rec[r] = (4*h[r]*local_R0*SSB[r]) / ((1-h[r])*S0[r] + (5*h[r]-1)*SSB[r]) # get local beverton holt
      }

      # Global Density Dependence
      if(recruitment_dd == 1) {
        rec[r] = (4*h[r]*R0*sum(SSB)) / ((1-h[r])*S0[r] + (5*h[r]-1)*sum(SSB)) * Rec_Prop[r] # get global beverton holt and then apportion to different regions
      }
    } # end r loop

  } # end Beverton-Holt

  return(rec)
}
