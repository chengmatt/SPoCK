#' Setup containers for simulation and output into global environment
#'
#' @param sim_list List object containing simulation components (should have n_regions, n_ages, n_sexes, n_sims, sim_list$n_fish_fleets, sim_list$n_srv_fleets)
#'
#' @export Setup_Sim_Containers
#'
Setup_Sim_Containers <- function(sim_list) {

  # Biological Containers
  sim_list$Init_NAA <- array(0, dim = c(sim_list$n_regions, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims))
  sim_list$Init_NAA_next_year <- sim_list$Init_NAA
  sim_list$NAA <- array(0, dim = c(sim_list$n_yrs+1, sim_list$n_regions, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims))
  sim_list$Z <- array(0, dim = c(sim_list$n_yrs, sim_list$n_regions, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims))
  sim_list$SSB <- array(0, dim = c(sim_list$n_yrs, sim_list$n_regions, sim_list$n_sims))
  sim_list$Total_Biom <- array(0, dim = c(sim_list$n_yrs, sim_list$n_regions, sim_list$n_sims))
  sim_list$ln_rec_devs <- array(0, dim = c(sim_list$n_yrs, sim_list$n_regions, sim_list$n_sims))

  # Fishery Containers
  sim_list$Obs_Catch <- array(0, dim = c(sim_list$n_yrs, sim_list$n_regions, sim_list$n_fish_fleets, sim_list$n_sims))
  sim_list$Obs_FishAgeComps <- array(0, dim = c(sim_list$n_yrs, sim_list$n_regions, sim_list$n_ages, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims))
  sim_list$CAA <- array(0, dim = c(sim_list$n_yrs, sim_list$n_regions, sim_list$n_ages, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims))
  sim_list$True_Catch <- array(0, dim = c(sim_list$n_yrs, sim_list$n_regions, sim_list$n_fish_fleets, sim_list$n_sims))
  sim_list$Obs_Catch <- array(0, dim = c(sim_list$n_yrs, sim_list$n_regions, sim_list$n_fish_fleets, sim_list$n_sims))

  # Survey Containers
  sim_list$Obs_SrvIdx <- array(0, dim = c(sim_list$n_yrs, sim_list$n_regions, sim_list$n_srv_fleets, sim_list$n_sims))
  sim_list$True_SrvIdx <- array(0, dim = c(sim_list$n_yrs, sim_list$n_regions, sim_list$n_srv_fleets, sim_list$n_sims))
  sim_list$Obs_SrvAgeComps <- array(0, dim = c(sim_list$n_yrs, sim_list$n_regions, sim_list$n_ages, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims))
  sim_list$Srv_IAA <- array(0, dim = c(sim_list$n_yrs, sim_list$n_regions, sim_list$n_ages, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims))

  return(sim_list)
}
