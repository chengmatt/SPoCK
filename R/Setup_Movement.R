#' Setup Movement Processes for SPoRC
#'
#' @param input_list List containing data, parameter, and map lists for the model.
#' @param do_recruits_move Integer flag (0 or 1) indicating whether recruits move.
#'   Default is 0 (do not move).
#' @param use_fixed_movement Integer flag (0 or 1) indicating whether to use a fixed movement matrix (1) or estimate movement parameters (0).
#'   Default is 0.
#' @param Fixed_Movement Numeric array for fixed movement matrix dimensioned by \code{[n_regions, n_regions, n_years, n_ages, n_sexes]}.
#'   Default is an array of ones.
#' @param Use_Movement_Prior Integer flag (0 or 1) indicating whether to use movement priors.
#'   Default is 0 (no priors).
#' @param Movement_prior Numeric vector or array specifying prior values for movement parameters.
#'   If a vector, a constant prior is applied across all dimensions.
#' @param Movement_ageblk_spec Either:
#'   \itemize{
#'     \item Character string \code{"constant"} for age-invariant movement, or
#'     \item A list of numeric vectors specifying age blocks sharing parameters.
#'   }
#'   For example, \code{list(c(1:6), c(7:10), c(11:n_ages))} defines three age blocks where:
#'   \itemize{
#'     \item ages 1 to 6 share parameters,
#'     \item ages 7 to 10 share parameters,
#'     \item ages 11 to \code{n_ages} share parameters.
#'   }
#'   To specify age-invariant movement, use either \code{"constant"} or \code{list(c(1:n_ages))}.
#' @param Movement_yearblk_spec Either:
#'   \itemize{
#'     \item Character string \code{"constant"} for time-invariant movement, or
#'     \item A list of numeric vectors specifying year blocks sharing movement parameters.
#'   }
#' @param Movement_sexblk_spec Either:
#'   \itemize{
#'     \item Character string \code{"constant"} for sex-invariant movement, or
#'     \item A list of numeric vectors specifying sex blocks sharing movement parameters.
#'   }
#' @param cont_vary_movement Character string specifying continuous varying movement type.
#'   Available options:
#'   \itemize{
#'     \item \code{"none"}
#'     \item \code{"iid_y"} (iid deviations by year)
#'     \item \code{"iid_a"} (iid deviations by age)
#'     \item \code{"iid_y_a"} (iid deviations by year and age)
#'     \item \code{"iid_y_s"} (iid deviations by year and sex)
#'     \item \code{"iid_a_s"} (iid deviations by age and sex)
#'     \item \code{"iid_y_a_s"} (iid deviations by year, age, and sex)
#'   }
#'   Default is \code{"none"}.
#' @param Movement_cont_pe_pars_spec Character string specifying process error parameter sharing.
#'   Available options:
#'   \itemize{
#'     \item \code{"est_shared_r"}
#'     \item \code{"est_shared_a"}
#'     \item \code{"est_shared_s"}
#'     \item \code{"est_shared_r_a"}
#'     \item \code{"est_shared_a_s"}
#'     \item \code{"est_shared_r_s"}
#'     \item \code{"est_shared_r_a_s"}
#'     \item \code{"est_all"}
#'     \item \code{"fix"}
#'     \item \code{"none"}
#'   }
#'   Default is \code{"none"}.
#' @param ... Additional parameters such as starting values for \code{move_pars}, \code{logit_move_devs}, and \code{move_pe_pars}.
#'
#' @export Setup_Mod_Movement
Setup_Mod_Movement <- function(input_list,
                               do_recruits_move = 0,
                               use_fixed_movement = 0,
                               Fixed_Movement = NA,
                               Use_Movement_Prior = 0,
                               Movement_prior = NULL,
                               Movement_ageblk_spec = NULL,
                               Movement_yearblk_spec = NULL,
                               Movement_sexblk_spec = NULL,
                               cont_vary_movement = 'none',
                               Movement_cont_pe_pars_spec = 'none',
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

  if(!cont_vary_movement %in% c("none", "iid_y", "iid_a", "iid_y_a", 'iid_y_s', 'iid_a_s', "iid_y_a_s")) stop('Options for continuous movement is not correctly specified. The options are none, iid_y, iid_a, iid_y_a, iid_y_s, iid_a_s, iid_y_a_s.')
  else collect_message("Continuous movement specification is: ", cont_vary_movement)

  if(!Movement_cont_pe_pars_spec %in% c('est_shared_r', 'est_shared_a', "est_shared_s",
                                        'est_shared_r_a', 'est_shared_a_s', "est_shared_r_s",
                                        'est_shared_r_a_s', 'est_all', 'fix', 'none')) stop('Options for continuous movement process error is not correctly specified. The options are est_shared_r, est_shared_a, est_shared_s, est_shared_r_a, est_shared_a_s, est_shared_r_s, est_shared_r_a_s, est_all, fix, none')
  else collect_message("Continuous movement process error specification is: ", Movement_cont_pe_pars_spec)

  # define compatbility movement rules
  compatibility_rules <- list(
    "none" = c("fix", "none"),  # no movement variation
    "iid_y" = c("fix", "est_shared_r_a_s"),  # must share ages and sexes
    "iid_a" = c("fix", "est_shared_r_a_s", "est_shared_r_s", "est_shared_a_s", "est_shared_s"),  # must share sexes
    "iid_y_s" = c("fix", "est_shared_r_a_s", "est_shared_a_s", "est_shared_r_a", "est_shared_a"),  # must share ages
    "iid_y_a" = c("fix", "est_shared_r_a_s", "est_shared_r_s", "est_shared_a_s", "est_shared_s"),  # must share sexes
    "iid_a_s" = c("fix", "est_shared_r", "est_shared_a", "est_shared_s", "est_shared_r_a", "est_shared_a_s", "est_shared_r_s", "est_shared_r_a_s"),  # all vary - no constraints
    "iid_y_a_s" = c("fix", "est_shared_r", "est_shared_a", "est_shared_s", "est_shared_r_a", "est_shared_a_s", "est_shared_r_s", "est_shared_r_a_s")  # all vary - no constraints
  )

  # Check compatibility
  if (!Movement_cont_pe_pars_spec %in% compatibility_rules[[cont_vary_movement]]) {
    compatible_specs <- compatibility_rules[[cont_vary_movement]]
    stop("Incompatible parameter combination. For cont_vary_movement = '", cont_vary_movement,
         "', Movement_cont_pe_pars_spec must be one of: ",
         paste(compatible_specs, collapse = ", "))
  }

  # Fixed effects blocks
  if(!is.null(Movement_ageblk_spec)) if(!typeof(Movement_ageblk_spec) %in% c("list", "character", NULL)) stop("Movement fixed effects age blocks are not correctly specified, it needs to be either a list object or set at 'constant'. For example, if we had 10 ages and wanted 2 age blocks, this would be list(c(1:5), c(6:10)) such that ages 1 - 5 are a block, and ages 6 - 10 are a block.")
  if(!is.null(Movement_yearblk_spec)) if(!typeof(Movement_yearblk_spec) %in% c("list", "character", NULL)) stop("Movement fixed effects year blocks are not correctly specified, it needs to be either a list object or set at 'constant'. For example, if we had 10 years and wanted 2 year blocks, this would be list(c(1:5), c(6:10)) such that years 1 - 5 are a block, and years 6 - 10 are a block.")
  if(!is.null(Movement_sexblk_spec)) if(!typeof(Movement_sexblk_spec) %in% c("list", "character", NULL)) stop("Movement fixed effects sex blocks are not correctly specified, it needs to be either a list object or set at 'constant'. For example, if we had 2 sexes and wanted sex-specific movement, this would be list(1, 2).")

  # define for continuous varying movement
  cont_move_map <- data.frame(type = c("none", "iid_y", "iid_a", "iid_y_a", 'iid_y_s', 'iid_a_s', "iid_y_a_s"), num = c(0:6))
  cont_vary_movement_val <- cont_move_map$num[cont_move_map$type == cont_vary_movement] # look for number corresponding to specified option

  # Input variables into data list
  input_list$data$do_recruits_move <- do_recruits_move
  input_list$data$use_fixed_movement <- use_fixed_movement
  input_list$data$Fixed_Movement <- Fixed_Movement
  input_list$data$Use_Movement_Prior <- Use_Movement_Prior
  Movement_prior_vals = ifelse(is.null(Movement_prior), rep(1, input_list$data$n_regions), Movement_prior) # set up movement prior values
  input_list$data$Movement_prior <- array(Movement_prior_vals, dim = c(input_list$data$n_regions, input_list$data$n_regions, length(input_list$data$years), length(input_list$data$ages), input_list$data$n_sexes))
  input_list$data$cont_vary_movement <- cont_vary_movement_val

  # Setup parameters list
  starting_values <- list(...) # get starting values if there are any

  # Movement Parameters
  if("move_pars" %in% names(starting_values)) input_list$par$move_pars <- starting_values$move_pars
  else input_list$par$move_pars <- array(0, dim = c(input_list$data$n_regions, input_list$data$n_regions - 1, length(input_list$data$years), length(input_list$data$ages), input_list$data$n_sexes))

  # Movement deviations
  if("logit_move_devs" %in% names(starting_values)) input_list$par$logit_move_devs <- starting_values$logit_move_devs
  else input_list$par$logit_move_devs <- array(0, c(input_list$data$n_regions, input_list$data$n_regions - 1, length(input_list$data$years), length(input_list$data$ages), input_list$data$n_sexes))

  # Movement process error parameters
  if("move_pe_pars" %in% names(starting_values)) input_list$par$move_pe_pars <- starting_values$move_pe_pars
  else input_list$par$move_pe_pars <- array(0, dim = c(input_list$data$n_regions, max(4, length(input_list$data$ages)), input_list$data$n_sexes)) # max 4 parameters or the ages

  # Setup mapping list
  map_Movement_Pars <- input_list$par$move_pars # initialize array with same dimensions as parameters
  map_logit_move_devs <- array(NA, dim = dim(input_list$par$logit_move_devs))
  map_move_pe_pars <- array(NA, dim = dim(input_list$par$move_pe_pars))

  # Setup dimensions
  n_regions_from <- dim(map_Movement_Pars)[1]
  n_regions_to <- dim(map_Movement_Pars)[2]

  # Whether or not recruits move
  age_start <- ifelse(do_recruits_move == 0 && length(input_list$data$ages) >= 2, 2, 1)

  # Set up mapping for fixed-effects of movement parameters (i.e., if blocking as fixed effects)
  # If movement is constant for either ages, years, or sexes
  if(is.character(Movement_ageblk_spec)){
    if(Movement_ageblk_spec == "constant") Movement_ageblk_spec_vals = list(input_list$data$ages)
  } else Movement_ageblk_spec_vals = Movement_ageblk_spec

  if(is.character(Movement_yearblk_spec)){
    if(Movement_yearblk_spec == "constant") Movement_yearblk_spec_vals = list(input_list$data$years)
  } else Movement_yearblk_spec_vals = Movement_yearblk_spec

  if(is.character(Movement_sexblk_spec)){
    if(Movement_sexblk_spec == "constant") Movement_sexblk_spec_vals = list(1:input_list$data$n_sexes)
  } else Movement_sexblk_spec_vals = Movement_sexblk_spec

  if(input_list$data$n_regions > 1) {
    counter <- 1
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

  # Set up mapping for logit move devs and process error parameters
  if(input_list$data$n_regions > 1 && input_list$data$cont_vary_movement > 0 && input_list$data$use_fixed_movement == 0) {
    counter <- 1

    for(r in 1:input_list$data$n_regions) {
      for(rr in 1:(input_list$data$n_regions - 1)) {

        if(cont_vary_movement %in% c('iid_y')) {
          for(y in 1:length(input_list$data$years)) {
            map_logit_move_devs[r,rr,y,,] <- counter
            counter <- counter + 1
          } # end y loop
        } # end if iid_y

        if(cont_vary_movement %in% c('iid_a')) {
          for(a in age_start:length(input_list$data$ages)) {
            map_logit_move_devs[r,rr,,a,] <- counter
            counter <- counter + 1
          } # end a loop
        } # end if iid_a

        if(cont_vary_movement %in% c('iid_y_a')) {
          for(y in 1:length(input_list$data$years)) {
            for(a in age_start:length(input_list$data$ages)) {
              map_logit_move_devs[r,rr,y,a,] <- counter
              counter <- counter + 1
            } # end a loop
          } # end y loop
        } # end if iid_y_a

        if(cont_vary_movement %in% c('iid_y_s')) {
          for(y in 1:length(input_list$data$years)) {
            for(s in 1:input_list$data$n_sexes) {
              map_logit_move_devs[r,rr,y,,s] <- counter
              counter <- counter + 1
            } # end s loop
          } # end y loop
        } # end if iid_y_s

        if(cont_vary_movement %in% c('iid_a_s')) {
          for(a in age_start:length(input_list$data$ages)) {
            for(s in 1:input_list$data$n_sexes) {
              map_logit_move_devs[r,rr,,a,s] <- counter
              counter <- counter + 1
            } # end s loop
          } # end y loop
        } # end if iid_a_s

        if(cont_vary_movement %in% c('iid_y_a_s')) {
          for(y in 1:length(input_list$data$years)) {
            for(a in age_start:length(input_list$data$ages)) {
              for(s in 1:input_list$data$n_sexes) {
                map_logit_move_devs[r,rr,y,a,s] <- counter
                counter <- counter + 1
              } # end s loop
            } # end a loop
          } # end y loop
        } # end if iid_y_a_s

      } # end rr
    } # end r loop

    # Mapping for movement process error deviations
    if(Movement_cont_pe_pars_spec %in% c("fix", "none")) map_move_pe_pars <- map_move_pe_pars
    if(Movement_cont_pe_pars_spec == 'est_all') map_move_pe_pars[] <- 1:length(map_move_pe_pars)
    if(Movement_cont_pe_pars_spec %in% c('est_shared_r', 'est_shared_a', "est_shared_s", 'est_shared_r_a', 'est_shared_a_s', 'est_shared_r_s', 'est_shared_r_a_s')) {
      counter <- 1
      for(r in 1:input_list$data$n_regions) {
        for(a in age_start:length(input_list$data$ages)) {
          for(s in 1:input_list$data$n_sexes) {
            # Sharing process error parameters across origin regions
            if(Movement_cont_pe_pars_spec == 'est_shared_r' && r == 1) {
              map_move_pe_pars[,a,s] <- counter
              counter <- counter + 1
            }
            # Sharing process error parameters across ages
            if(Movement_cont_pe_pars_spec == 'est_shared_a' && a == 1) {
              map_move_pe_pars[r,,s] <- counter
              counter <- counter + 1
            }
            # Sharing process error parameters across sexes
            if(Movement_cont_pe_pars_spec == 'est_shared_s' && s == 1) {
              map_move_pe_pars[r,a,] <- counter
              counter <- counter + 1
            }
            # Sharing process error parameters across regions and ages
            if(Movement_cont_pe_pars_spec == 'est_shared_r_a' && r == 1 && a == 1) {
              map_move_pe_pars[,,s] <- counter
              counter <- counter + 1
            }
            # Sharing process error parameters across ages and sexes
            if(Movement_cont_pe_pars_spec == 'est_shared_a_s' && a == 1 && s == 1) {
              map_move_pe_pars[r,,] <- counter
              counter <- counter + 1
            }
            # Sharing process error parameters across regions and sexes
            if(Movement_cont_pe_pars_spec == 'est_shared_r_s' && r == 1 && s == 1) {
              map_move_pe_pars[,a,] <- counter
              counter <- counter + 1
            }
            # Sharing process error parameters across regions, ages, and sexes
            if(Movement_cont_pe_pars_spec == 'est_shared_r_a_s' && r == 1 && a == 1 && s == 1) {
              map_move_pe_pars[,,] <- counter
              counter <- counter + 1
            }
          } # end s loop
        } # end a loop
      } # end r loop
    }
  }  else {
    map_logit_move_devs <- factor(rep(NA, length(input_list$par$logit_move_devs)))
    map_move_pe_pars <- factor(rep(NA, length(input_list$par$move_pe_pars)))
  }

  # Print out movement specifications and error handling
  if(is.list(Movement_sexblk_spec)) collect_message("Movement fixed effect blocks are specified with ", length(Movement_sexblk_spec), " sex blocks") else collect_message("Movement fixed effect blocks are sex-invariant")
  if(is.list(Movement_yearblk_spec)) collect_message("Movement fixed effect blocks are specified with ", length(Movement_yearblk_spec), " year blocks") else collect_message("Movement fixed effect blocks are time-invariant")
  if(is.list(Movement_ageblk_spec)) collect_message("Movement fixed effect blocks are specified with ", length(Movement_ageblk_spec), " age blocks") else collect_message("Movement fixed effect blocks are age-invariant")

  # Input into data and mapping list
  input_list$map$move_pars <- factor(map_Movement_Pars)
  input_list$map$logit_move_devs <- factor(map_logit_move_devs)
  input_list$map$move_pe_pars <- factor(map_move_pe_pars)
  input_list$data$map_Movement_Pars <- array(as.numeric(input_list$map$move_pars), dim = dim(input_list$par$move_pars))
  input_list$data$map_logit_move_devs <- array(as.numeric(input_list$map$logit_move_devs), dim = dim(input_list$par$logit_move_devs))

  # Print all messages if verbose is TRUE
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}
