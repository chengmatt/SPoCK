#' Get Deterministic Recruitment
#'
#' @param recruitment_model == 0, mean recruitment, == 1 beverton holt recruitment with steepness
#' @param R0 virgin or mean recruitment
#' @param h steepness
#' @param n_ages number of ages
#' @param WAA vector of weight at age
#' @param MatAA vector of maturity at age
#' @param natmort vector of natural mortality at age
#' @param SSB_vals spawning stock biomass value
#' @param y year
#' @param rec_lag recruitment lag for indexing SSB year
#'
#' @returns Determinstic recruitment values from mean recruitment or beverton holt
#' @keywords internal
#'
Get_Det_Recruitment <- function(recruitment_model,
                                y,
                                rec_lag,
                                R0,
                                h,
                                n_ages,
                                WAA,
                                MatAA,
                                natmort,
                                SSB_vals
                                ) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  if(recruitment_model == 0) rec = R0 # mean recruitment
  if(recruitment_model == 1) { # beverton-holt
    # Calculate unexploited naa per recruit
    tmp_naa = rep(0, n_ages)
    tmp_naa[1] = 1
    for(a in 2:n_ages) tmp_naa[a] = tmp_naa[a-1] * exp(-natmort[a-1])
    S0 = sum(tmp_naa * WAA * MatAA) * R0 # get unexploited ssb
    if(y <= rec_lag) SSB = S0 else SSB = SSB_vals[y-rec_lag] # use SSB indexed by recruitment lag year
    rec = (4*R0*h*SSB) / ((1-h)*R0*(S0/R0) + (5*h - 1)*SSB)
  } # end if for beveteron holt
  return(rec)
}
