#' Setup containers for simulation and output into global environment
#'
#' @param n_sims Number of simulations
#' @param n_yrs Number of years
#' @param n_regions Number of regions
#' @param n_ages Number of ages
#' @param n_sexes Number of sexes
#' @param n_fish_fleets Number of fishery fleets
#' @param n_srv_fleets Number of survey fleets
#'
#' @export Setup_Sim_Containers
#'
Setup_Sim_Containers <- function(n_sims = n_sims,
                                 n_yrs = n_yrs,
                                 n_regions = n_regions,
                                 n_ages = n_ages,
                                 n_sexes = n_sexes,
                                 n_fish_fleets = n_fish_fleets,
                                 n_srv_fleets = n_srv_fleets) {

  # Biological Containers
  Init_NAA <<- array(0, dim = c(n_regions, n_ages, n_sexes, n_sims))
  Init_NAA_next_year <<- Init_NAA
  NAA <<- array(0, dim = c(n_yrs+1, n_regions, n_ages, n_sexes, n_sims))
  Z <<- array(0, dim = c(n_yrs, n_regions, n_ages, n_sexes, n_sims))
  SSB <<- array(0, dim = c(n_yrs, n_regions, n_sims))
  Total_Biom <<- array(0, dim = c(n_yrs, n_regions, n_sims))
  ln_rec_devs <<- array(0, dim = c(n_yrs, n_regions, n_sims))

  # Fishery Containers
  Obs_Catch <<- array(0, dim = c(n_yrs, n_regions, n_fish_fleets, n_sims))
  Obs_FishAgeComps <<- array(0, dim = c(n_yrs, n_regions, n_ages, n_sexes, n_fish_fleets, n_sims))
  CAA <<- array(0, dim = c(n_yrs, n_regions, n_ages, n_sexes, n_fish_fleets, n_sims))
  True_Catch <<- array(0, dim = c(n_yrs, n_regions, n_fish_fleets, n_sims))
  Obs_Catch <<- array(0, dim = c(n_yrs, n_regions, n_fish_fleets, n_sims))

  # Survey Containers
  Obs_SrvIdx <<- array(0, dim = c(n_yrs, n_regions, n_srv_fleets, n_sims))
  True_SrvIdx <<- array(0, dim = c(n_yrs, n_regions, n_srv_fleets, n_sims))
  Obs_SrvAgeComps <<- array(0, dim = c(n_yrs, n_regions, n_ages, n_sexes, n_srv_fleets, n_sims))
  Srv_IAA <<- array(0, dim = c(n_yrs, n_regions, n_ages, n_sexes, n_srv_fleets, n_sims))

}
