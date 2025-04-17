#' Set up simulation dimensions
#'
#' @param n_sims Number of simulations
#' @param n_yrs Number of years
#' @param n_regions Number of regions
#' @param n_ages Number of ages
#' @param n_sexes Number of sexes
#' @param n_fish_fleets Number of fishery fleets
#' @param n_srv_fleets Number of survey fleets
#'
#' @export Setup_Sim_Dim
#'
Setup_Sim_Dim <- function(n_sims,
                          n_yrs,
                          n_regions,
                          n_ages,
                          n_sexes,
                          n_fish_fleets,
                          n_srv_fleets) {

  # ouput variables into global environment
  n_sims <<- n_sims
  n_yrs <<- n_yrs
  n_regions <<- n_regions
  n_ages <<- n_ages
  n_sexes <<- n_sexes
  n_fish_fleets <<- n_fish_fleets
  n_srv_fleets <<- n_srv_fleets
  init_iter <<- n_ages * 10

}

#' Set up model dimensions
#'
#' @param n_regions Number of regions
#' @param ages vector of ages
#' @param n_sexes Number of sexes
#' @param n_fish_fleets Number of fishery fleets
#' @param n_srv_fleets Number of survey fleets
#' @param years vector of years
#' @param lens vector of lengths (can just input 1 if not fitting lengths)
#'
#' @returns Returns a list object of length 3, with a data list, a par list, and a map list
#' @export Setup_Mod_Dim
#'
Setup_Mod_Dim <- function(years,
                          ages,
                          lens,
                          n_regions,
                          n_sexes,
                          n_fish_fleets,
                          n_srv_fleets) {

  # Create empty list
  input_list <- list(data = list(), par = list(), map = list())

  # ouput variables into list
  input_list$data$years <- years
  input_list$data$n_regions <- n_regions
  input_list$data$ages <- ages
  input_list$data$lens <- lens
  input_list$data$n_sexes <- n_sexes
  input_list$data$n_fish_fleets <- n_fish_fleets
  input_list$data$n_srv_fleets <- n_srv_fleets

  message("Number of Years: ", length(years))
  message("Number of Regions: ", n_regions)
  message("Number of Age Bins: ", length(ages))
  message("Number of Length Bins: ", length(lens))
  message("Number of Sexes: ", n_sexes)
  message("Number of Fishery Fleets: ", n_fish_fleets)
  message("Number of Survey Fleets: ", n_srv_fleets)

  return(input_list)

}
