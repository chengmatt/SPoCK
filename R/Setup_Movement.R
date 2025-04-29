#' Setup model movement processes
#'
#' @param input_list List containing a data list, parameter list, and map list
#' @param do_recruits_move Whether recruits move, 1 (they move), 0 (they don't move), default is 0
#' @param use_fixed_movement Whether to use a fixed movement matrix, == 0 don't use, == 1 use
#' @param Fixed_Movement Input for a fixed movement matrix dimensioned by n_regions, n_regions, n_years, n_ages, n_sexes
#' @param Use_Movement_Prior Whether or not to use movement priors == 0 don't use, == 1 use
#' @param Movement_prior Movement prior values, can be dimension as a vector of (results in a constant prior across all regions from, years, ages, and sexes), or as an array by n_regions, n_regions, n_years, n_ages, n_sexes
#' @param ... Additional inputs for starting values for movement parameters (move_pars)
#' @param Movement_ageblk_spec Either a character string specifiying "constant" age movement, or a ist object specifying which ages to group together and block (i.e., which ages to share parameters for). Each element in the list should be a vector specifying the range of ages for which to share or block parameters for. For example, list(c(1:6), c(7:10), c(11:n_ages)) has 3 age blocks, where ages 1-6 have the same paraemters, 7:10 have the same parameters (but diffeerent from 1:6), and 11:n_ages haave the same parameters (but different from 1:6 and 7:10).
#' The default behavior is to estimate unique movement parmaeters for all ages. To estimate age-invariant movement this list would be specified as: list(c(1:n_ages)).
#' @param Movement_yearblk_spec Either a character string specifiying "constant" year movement, list object specifying which years to group together and block (i.e., which years to share parameters for). Each element in the list should be a vector specifying the range of years for which to share or block parameters for. For example, list(c(1:6), c(7:10), c(11:n_years)) has 3 year blocks, where years 1-6 have the same paraemters, 7:10 have the same parameters (but diffeerent from 1:6), and 11:n_years haave the same parameters (but different from 1:6 and 7:10).
#' The default behavior is to estimate unique movement parameters for all years To estimate time-invariant movement this list would be specified as: list(c(1:n_years)).
#' @param Movement_sexblk_spec Either a character string specifiying "constant" sex movement, list object specifying which sexes to group together and block (i.e., which sexes to share parameters for). Each element in the list should be a vector specifying the range of sexes for which to share or block parameters for. For example, list(1, 2) has 2 sex blocks, where each sex has unique parameters
#' The default behavior is to estimate unique movement parameters for all sexes. To estimate sex-invariant movement this list would be specified as: list(c(1:n_sexes)).
#'
#' @export Setup_Mod_Movement
#'
Setup_Mod_Movement <- function(input_list,
                               do_recruits_move = 0,
                               use_fixed_movement = 0,
                               Fixed_Movement = NA,
                               Use_Movement_Prior = 0,
                               Movement_prior = NULL,
                               Movement_ageblk_spec = NULL,
                               Movement_yearblk_spec = NULL,
                               Movement_sexblk_spec = NULL,
                               ...
                               ) {

  messages_list <<- character(0) # string to attach to for printing messages

  # if no fixed movement matrix provided
  if(is.na(sum(Fixed_Movement))) {
    Fixed_Movement <- array(1, dim = c(input_list$data$n_regions, input_list$data$n_regions, length(input_list$data$years), length(input_list$data$ages), input_list$data$n_sexes))
  }

  # Check fixed movement matrix
  if(use_fixed_movement == 1) check_data_dimensions(Fixed_Movement, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, what = 'Fixed_Movement')

  if(!do_recruits_move %in% c(0,1)) stop('Movement for recruits is not correctly specified. The options are do_recruits_move == 0 (they dont move), or == 1 (they move)')
  else collect_message("Recruits are: ", ifelse(do_recruits_move == 0, "Not Moving", "Moving"))

  if(!use_fixed_movement %in% c(0,1)) stop('Options for fixing movement are not correctly specified. The options are use_fixed_movement == 0 (dont use and estiamte movement parameters), or == 1 (use)')
  else collect_message("Movement is: ", ifelse(use_fixed_movement == 0, "Estimated", "Fixed"))

  if(!Use_Movement_Prior %in% c(0,1)) stop('Options for movement priors not correctly specified. The options are Use_Movement_Prior == 0 (dont use), or == 1 (use)')
  else collect_message("Movement priors are: ", ifelse(Use_Movement_Prior == 0, "Not Used", "Used"))

  if(!is.null(Movement_ageblk_spec)) if(!typeof(Movement_ageblk_spec) %in% c("list", "character", NULL)) stop("Movement age blocks are not correctly specified, it needs to be either a list object or set at 'constant'. For example, if we had 10 ages and wanted 2 age blocks, this would be list(c(1:5), c(6:10)) such that ages 1 - 5 are a block, and ages 6 - 10 are a block.")
  if(!is.null(Movement_yearblk_spec)) if(!typeof(Movement_yearblk_spec) %in% c("list", "character", NULL)) stop("Movement year blocks are not correctly specified, it needs to be either a list object or set at 'constant'. For example, if we had 10 years and wanted 2 year blocks, this would be list(c(1:5), c(6:10)) such that years 1 - 5 are a block, and years 6 - 10 are a block.")
  if(!is.null(Movement_sexblk_spec)) if(!typeof(Movement_sexblk_spec) %in% c("list", "character", NULL)) stop("Movement sex blocks are not correctly specified, it needs to be either a list object or set at 'constant'. For example, if we had 2 sexes and wanted sex-specific movement, this would be list(1, 2).")

  # Input variables into data list
  input_list$data$do_recruits_move <- do_recruits_move
  input_list$data$use_fixed_movement <- use_fixed_movement
  input_list$data$Fixed_Movement <- Fixed_Movement
  input_list$data$Use_Movement_Prior <- Use_Movement_Prior
  Movement_prior_vals = ifelse(is.null(Movement_prior), rep(1, input_list$data$n_regions), Movement_prior) # set up movement prior values
  input_list$data$Movement_prior <- array(Movement_prior_vals, dim = c(input_list$data$n_regions, input_list$data$n_regions, length(input_list$data$years), length(input_list$data$ages), input_list$data$n_sexes))

  # Setup parameters list
  starting_values <- list(...) # get starting values if there are any

  # Movement Parameters
  if("move_pars" %in% names(starting_values)) input_list$par$move_pars <- starting_values$move_pars
  else input_list$par$move_pars <- array(0, dim = c(input_list$data$n_regions, input_list$data$n_regions - 1, length(input_list$data$years), length(input_list$data$ages), input_list$data$n_sexes))

  # Setup mapping list
  map_Movement_Pars <- input_list$par$move_pars # initialize array with same dimensions as parameters

  # Setup dimensions
  n_regions_from <- dim(map_Movement_Pars)[1]
  n_regions_to <- dim(map_Movement_Pars)[2]

  # Initialize counter to loop through (specify map as NA if using fixed movement)
  if(use_fixed_movement == 0) counter <- 1
  if(use_fixed_movement == 1 || input_list$data$n_regions == 1) counter <- NA

  # If movement is constant for either ages, years, or sexes
  if(is.character(Movement_ageblk_spec)){
    if(Movement_ageblk_spec == "constant") Movement_ageblk_spec_vals = list(input_list$data$ages)
  } else Movement_ageblk_spec_vals = Movement_ageblk_spec

  if(is.character(Movement_yearblk_spec)){
    if(Movement_yearblk_spec == "constant") Movement_yearblk_spec_vals = list(input_list$data$years)
  } else Movement_yearblk_spec_vals = Movement_yearblk_spec

  if(is.character(Movement_sexblk_spec)){
    if(Movement_sexblk_spec == "constant") Movement_sexblk_spec_vals = list(1:input_list$data$n_sexes)
  } else Movement_sexblk_spec_vals = Movement_yearblk_spec

  counter <- 1
  if(input_list$data$n_regions > 1) {
    for(ageblk in 1:length(Movement_ageblk_spec_vals)) {
      map_a <- Movement_ageblk_spec_vals[[ageblk]] # get ages to block and map off
      for(yearblk in 1:length(Movement_yearblk_spec_vals)) {
        map_y <- Movement_yearblk_spec_vals[[yearblk]] # get years to block and map off
        for(sexblk in 1:length(Movement_sexblk_spec_vals)) {
          map_s <- Movement_sexblk_spec_vals[[sexblk]] # get sexes to block and map off

          # Now, loop through each combination and increment get unique indices
          map_idx <- array(0, dim = c(n_regions_from, n_regions_to))
          for(i in 1:n_regions_from) {
            for(j in 1:n_regions_to) {
              map_idx[i,j] <- counter
              counter <- counter + 1 # increment counter
            } # end j loop
          } # end i loop

          # Now input these unique indices into the map array
          for(a in map_a) for(y in map_y) for(s in map_s) map_Movement_Pars[,,y,a,s] <- map_idx

        } # end sex block
      } # end year block
    } # end age block
  } else map_Movement_Pars <- factor(rep(NA, length(input_list$par$move_pars)))

  # Print out movement specifications and error handling
  if(is.list(Movement_sexblk_spec)) collect_message("Movement is specified with ", length(Movement_sexblk_spec), " sex blocks") else collect_message("Movement is sex-invariant")
  if(is.list(Movement_yearblk_spec)) collect_message("Movement is specified with ", length(Movement_yearblk_spec), " year blocks") else collect_message("Movement is time-invariant")
  if(is.list(Movement_ageblk_spec)) collect_message("Movement is specified with ", length(Movement_ageblk_spec), " age blocks") else collect_message("Movement is age-invariant")

  # Input into data and mapping list
  input_list$map$move_pars <- factor(map_Movement_Pars)
  input_list$data$map_Movement_Pars <- array(as.numeric(input_list$map$move_pars), dim = dim(input_list$par$move_pars))

  # Print all messages if verbose is TRUE
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}
