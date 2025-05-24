#' Set up simulation dimensions
#'
#' @param n_sims Number of simulations
#' @param n_yrs Number of years
#' @param n_regions Number of regions
#' @param n_ages Number of ages
#' @param n_sexes Number of sexes
#' @param n_fish_fleets Number of fishery fleets
#' @param n_srv_fleets Number of survey fleets
#' @param run_feedback Boolean for whether to run or not run feedback management loop
#' @param feedback_start_yr If MSE is run, when is the first year feedback starts
#'
#' @export Setup_Sim_Dim
#'
Setup_Sim_Dim <- function(n_sims,
                          n_yrs,
                          n_regions,
                          n_ages,
                          n_sexes,
                          n_fish_fleets,
                          n_srv_fleets,
                          run_feedback,
                          feedback_start_yr
                          ) {

  sim_list <- list() # setup empty list

  # output dimensions into list
  sim_list$n_sims <- n_sims
  sim_list$n_yrs <- n_yrs
  sim_list$n_regions <- n_regions
  sim_list$n_ages <- n_ages
  sim_list$n_sexes <- n_sexes
  sim_list$n_fish_fleets <- n_fish_fleets
  sim_list$n_srv_fleets <- n_srv_fleets
  sim_list$init_iter <- n_ages * 10
  sim_list$feedback_start_yr <- feedback_start_yr
  sim_list$run_feedback <- run_feedback

  return(sim_list)

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
#' @param verbose Whether to print messages
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
                          n_srv_fleets,
                          verbose = FALSE
                          ) {

  messages_list <<- character(0) # string to attach to for printing messages

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
  input_list$verbose <- verbose

  collect_message("Number of Years: ", length(years))
  collect_message("Number of Regions: ", n_regions)
  collect_message("Number of Age Bins: ", length(ages))
  collect_message("Number of Length Bins: ", length(lens))
  collect_message("Number of Sexes: ", n_sexes)
  collect_message("Number of Fishery Fleets: ", n_fish_fleets)
  collect_message("Number of Survey Fleets: ", n_srv_fleets)

  # Print all messages if verbose is TRUE
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)

}
