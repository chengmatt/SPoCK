#' Set up recruitment dynamics for simulation
#'
#' @param do_recruits_move whether recruits move. Character string either specified as 'dont_move', or 'move'
#' @param base_rec_sexratio base recruitment sex-ratio value
#' @param rec_sexratio_vary whether recruitment sex-ratio varies. Options include: constant
#' @param base_r0 base R0 or mean recruitment value
#' @param r0_vary whether r0 or mean recruitment varies. Options include: constant
#' @param base_h base steepness value
#' @param init_sigmaR Sigma R for initial devs
#' @param sigmaR Sigma R for everything else
#' @param recruitment_opt character string as either "mean_rec" or "bh_rec"
#' @param rec_dd recruitment density dependence; character string as either "global" or "local"
#' @param sim_list Simulation list
#' @param init_dd initial age density dependence; character string as either "global" or "local"
#' @param rec_lag Numeric, recruitment lag value
#'
#' @export Setup_Sim_Rec
#'
Setup_Sim_Rec <- function(
    do_recruits_move,
    base_rec_sexratio,
    rec_sexratio_vary,
    base_r0,
    r0_vary,
    base_h,
    init_sigmaR,
    sigmaR,
    recruitment_opt,
    rec_dd,
    init_dd,
    sim_list,
    rec_lag
) {


  if(do_recruits_move == 'dont_move') sim_list$do_recruits_move <- 0
  if(do_recruits_move == 'move') sim_list$do_recruits_move <- 1
  if(sim_list$do_recruits_move == 0) sim_list$move_age <- 2 else sim_list$move_age <- 1 # what age to start movement of individuals

  if(recruitment_opt == "mean_rec") sim_list$recruitment_opt <- 0
  if(recruitment_opt == "bh_rec") sim_list$recruitment_opt <- 1

  if(rec_dd == "global") sim_list$rec_dd <- 0
  if(rec_dd == "local") sim_list$rec_dd <- 1
  if(init_dd == "global") sim_list$init_dd <- 0
  if(init_dd == "local") sim_list$init_dd <- 1

  # recruitment variability
  sim_list$init_sigmaR <- array(init_sigmaR, dim = c(sim_list$n_regions, 1))
  sim_list$sigmaR <- array(sigmaR, dim = c(sim_list$n_regions, sim_list$n_yrs))

  # Set up containers for r0, h, init_sigmaR, sigmaR, and recruitment sex ratio
  r0 <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims))
  h <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims))
  rec_sexratio <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_sexes, sim_list$n_sims))

  for(sim in 1:sim_list$n_sims) {
    for(r in 1:sim_list$n_regions) {
      for(y in 1:sim_list$n_yrs) {

        # Fill in values
        h[r,y,sim] <- base_h[r] # fill in steepness
        if(r0_vary == "constant") r0[r,y,sim] <- base_r0[r] # fill in r0 constant
        if(rec_sexratio_vary == "constant") for(s in 1:sim_list$n_sexes) rec_sexratio[r,y,s,sim] <- base_rec_sexratio[s] # fill in constant recruitment sex-ratio

      } # end y loop
    } # end r loop
  } # end sim loop

  # output these into environment
  sim_list$h <- h
  sim_list$r0 <- r0
  sim_list$rec_sexratio <- rec_sexratio
  sim_list$rec_lag <- rec_lag

  return(sim_list)

}



#' Setup model objects for specifying recruitment module and associated processes
#'
#' @param input_list List containing a data list, parameter list, and map list
#' @param rec_model recruitment model, options are mean_rec, bh_rec
#' @param rec_lag duration of recruitment lag with SSB
#' @param Use_h_prior whether or not to use a steepness prior == 0 (don't use), == 1 (use)
#' @param h_mu vector of steepness mean if priors are used (dimensioned by n_region)
#' @param h_sd vector of steepness sd if priors are used (dimensioned by n_region)
#' @param do_rec_bias_ramp whether or not to do recruitment bias ramp == 0 (don't use), == 1(use)
#' @param bias_year vector of years in which the bias ramp is invoked and when it changes
#' @param sigmaR_switch year in which sigmaR switches from an early value to late (set to 0 if no switching occurs)
#' @param sexratio vector of recruitment sex ratio dimensioned by n_sexes
#' @param init_age_strc whether initial age structure is initialized via == 0 (iteration) or via == 1 (geometric series)
#' @param init_F_prop Initial F proportion relative to a mean F value for initializing age structure
#' @param ... Additional arguments for starting values for recruitment parameters (ln_global_R0, Rec_prop, h, ln_InitDevs, ln_RecDevs, ln_sigmaR)
#' @param sigmaR_spec Specification for sigmaR. Default is NULL such that it is estimated for both early and late periods. Other options include "est_shared" which estiamtes sigmaR but shares the same value between early and late period and "fix" which fixes both values.
#' @param InitDevs_spec Specification for initial age dviations. Default is NULL such that is its estimated for all ages and regions. Other options include "est_shared_r" which estimates deviations for all ages but shares them across regions, or "fix" which fixeds all values.
#' @param RecDevs_spec Specification for recruitment deviations. Default is NULL such that it is estimated for all regions and years. Other options include "est_shared_r" which estimates deviations for all years, but shares them across regions (global recruitment deviaitons), or "fix" which fixes all values.
#' @param dont_est_recdev_last Numeric indicating the last x years to not estimate rec devs for. Default is 0.
#' @param h_spec Specification for steepness. Default is NULL, such that it is estimated for all reigons if rec_model == 1 (Beverton-Holt). Other options include "est_shared_r" which estimates h but shares across regions, or "fix" which fixes all steepness values. If rec_model == 0, h is fixed.
#' @param rec_dd Recruitment density dependence, Options are "local", "global", or NULL
#' @param t_spawn Spawn timing fraction
#'
#' @export Setup_Mod_Rec
#'
Setup_Mod_Rec <- function(input_list,
                          rec_model,
                          rec_dd = NULL,
                          rec_lag = 1,
                          Use_h_prior = 0,
                          h_mu = NA,
                          h_sd = NA,
                          do_rec_bias_ramp = 0,
                          bias_year = NA,
                          sigmaR_switch = 1,
                          dont_est_recdev_last = 0,
                          sexratio = 1,
                          init_age_strc = 0,
                          init_F_prop = 0,
                          sigmaR_spec = NULL,
                          InitDevs_spec = NULL,
                          RecDevs_spec = NULL,
                          h_spec = NULL,
                          t_spawn = 0,
                          ...
) {

  messages_list <<- character(0) # string to attach to for printing messages

  # Specify recruitment options
  if(rec_model == "mean_rec") {
    rec_model_val <- 0
    rec_dd_val <- 999
  }

  if(rec_model == "bh_rec") {
    rec_model_val <- 1
    if(is.null(rec_dd)) {
      rec_dd_val <- 1
    } else {
      if(rec_dd == 'local') rec_dd_val <- 0
      if(rec_dd == 'global') rec_dd_val <- 1
    }
  }

  if(!rec_model %in% c("mean_rec", "bh_rec")) stop("Please specify a valid recruitment form. These include: mean_rec or bh_rec")
  collect_message("Recruitment is specified as: ", rec_model)

  if(!is.null(rec_dd)) if(!rec_dd %in% c("local", "global")) stop("Please specify a valid recruitment density dependence form. These include: local or global")
  collect_message("Recruitment Density Dependence is specified as: ", rec_dd)

  if(rec_model != 'mean_rec') collect_message("Recruitment and SSB lag is specified as: ", rec_lag)
  if(rec_model == 'bh_rec') if(!Use_h_prior %in% c(0,1)) stop("Steepness priors are not specified as either: 0 (don't turn on), or 1 (turn on)")
  else collect_message("Steepness priors are: ", ifelse(Use_h_prior == 0, 'Not Used', 'Used'))

  if(!do_rec_bias_ramp %in% c(0,1)) stop("Recruitment bias ramp options are either: 0 (don't turn on), or 1 (turn on)")
  collect_message("Recruitment Bias Ramp is: ", ifelse(do_rec_bias_ramp == 0, "Off", 'On'))

  if(!init_age_strc %in% c(0,1)) stop("Initial Age Structure options are either: 0 (iterate), or 1 (geometric series solution)")
  collect_message("Initial Age Structure is: ", ifelse(init_age_strc == 0, "Iterated", 'Geometric Series Solution'))

  if(!is.numeric(sigmaR_switch)) stop("sigmaR_switch needs to be a numeric value. If you want to use a single sigmaR, please specify at 1")
  if(sigmaR_switch > 1) collect_message("Sigma R switches from an early period value to a late period value at year: ", sigmaR_switch)

  if(dont_est_recdev_last == 0) collect_message("Recruitment deviations for every year is estiamted")
  else collect_message("Recruitment deviations are not estimated for terminal year - ", dont_est_recdev_last, ". Recruitment during those periods are specified as the mean / deterministic recruitment")

  # Set up parameter list
  starting_values <- list(...) # get starting values if there are any

  # input variables into data list
  input_list$data$rec_model <- rec_model_val
  input_list$data$rec_dd <- rec_dd_val
  input_list$data$rec_lag <- rec_lag
  input_list$data$Use_h_prior <- Use_h_prior
  input_list$data$h_mu <- h_mu
  input_list$data$h_sd <- h_sd
  input_list$data$do_rec_bias_ramp <- do_rec_bias_ramp
  input_list$data$bias_year <- bias_year
  input_list$data$sigmaR_switch <- sigmaR_switch
  input_list$data$sexratio <- sexratio
  input_list$data$init_age_strc <- init_age_strc
  input_list$data$init_F_prop <- init_F_prop
  input_list$data$t_spawn <- t_spawn

  # Global R0
  if("ln_global_R0" %in% names(starting_values)) input_list$par$ln_global_R0 <- starting_values$ln_global_R0
  else input_list$par$ln_global_R0 <- log(15)

  # R0 proportion
  if("Rec_prop" %in% names(starting_values)) input_list$par$Rec_prop <- starting_values$Rec_prop
  else input_list$par$Rec_prop <- array(rep(0, input_list$data$n_regions - 1))

  # Steepness in bounded logit space (0.2 and 1)
  if("steepness_h" %in% names(starting_values)) input_list$par$steepness_h <- starting_values$steepness_h
  else input_list$par$steepness_h <- rep(0, input_list$data$n_regions)

  # Initial age deviations
  if("ln_InitDevs" %in% names(starting_values)) input_list$par$ln_InitDevs <- starting_values$ln_InitDevs
  else input_list$par$ln_InitDevs <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$ages) - 2))

  # Recruitment deviations
  if("ln_RecDevs" %in% names(starting_values)) input_list$par$ln_RecDevs <- starting_values$ln_RecDevs
  else input_list$par$ln_RecDevs <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years) - dont_est_recdev_last))

  # Recruitment variability (early period 1st element, late period 2nd element)
  if("ln_sigmaR" %in% names(starting_values)) input_list$par$ln_sigmaR <- starting_values$ln_sigmaR
  else input_list$par$ln_sigmaR <- c(0,0)

  # Setup mapping stuff
  # Standard deviation for initial age deviations and recruitment deviations
  if(!is.null(sigmaR_spec)) {
    # Share early and late period sigmaR and estimate
    if(sigmaR_spec == "est_shared") input_list$map$ln_sigmaR <- factor(c(1,1))
    if(sigmaR_spec == "fix_early_est_late") input_list$map$ln_sigmaR <- factor(c(NA, 1))
    # Fix both sigmaRs at starting values
    if(sigmaR_spec == "fix") input_list$map$ln_sigmaR <- factor(c(NA, NA))
    if(!sigmaR_spec %in% c("est_shared", "fix_early_est_late", "fix"))  stop("Please specify a valid recruitment variability option. These include: fix, fix_early_est_late, est_shared. Conversely, leave at NULL to estimate all recruitment variability parameters (early and late period).")
    else collect_message("Recruitment Variability is specified as: ", sigmaR_spec)
  } else collect_message("Recruitment Variability is estimated for both early and late periods")

  # Initial age deviations
  if(!is.null(InitDevs_spec)) {
    map_InitDevs <- input_list$par$ln_InitDevs # set up mapping for initial age deviations
    if(rec_dd == 'global' && InitDevs_spec != "est_shared_r" && input_list$data$n_regions > 1) stop("Please specify a valid initial age deviations option for global recruitment density dependence (should be est_shared_r or leave as NULL)!")
    # Share across regions and estimate
    if(InitDevs_spec == "est_shared_r") {
      for(r in 1:input_list$data$n_regions) map_InitDevs[r,] <- 1:length(map_InitDevs[1,]) # share parameters across regions
      input_list$map$ln_InitDevs <- factor(map_InitDevs)
    } # end if
    # Fix all initial deviations
    if(InitDevs_spec == "fix") input_list$map$ln_InitDevs <- factor(rep(NA, prod(dim(map_InitDevs))))
    if(!InitDevs_spec %in% c("est_shared_r", "fix"))  stop("Please specify a valid initial deviations option. These include: fix, est_shared_r. Conversely, leave at NULL to estimate all initial deviations.")
    else collect_message("Initial Deviations is specified as: ", InitDevs_spec)
  } else collect_message("Initial Age Deviations is estimated for all dimensions")

  # Recruitment deviations
  if(!is.null(RecDevs_spec)) {
    map_RecDevs <- input_list$par$ln_RecDevs # set up mapping for recruitment deviations
    if(rec_dd == 'global' && RecDevs_spec != "est_shared_r" && input_list$data$n_regions > 1) stop("Please specify a valid recruitment deviations option for global recruitment density dependence (should be est_shared_r)!")
    # Share across regions and estimate
    if(RecDevs_spec == "est_shared_r") {
      for(r in 1:input_list$data$n_regions) map_RecDevs[r,] <- 1:length(map_RecDevs[1,]) # share parameters across regions
      input_list$map$ln_RecDevs <- factor(map_RecDevs)
    } # end if
    # Fix all recruitment deviations
    if(RecDevs_spec == "fix") input_list$map$ln_RecDevs <- factor(rep(NA, prod(dim(map_RecDevs))))
    if(!RecDevs_spec %in% c("est_shared_r", "fix"))  stop("Please specify a valid recruitment deviations option. These include: fix, est_shared_r. Conversely, leave at NULL to estimate all recruitment deviations.")
    else collect_message("Recruitment Deviations is specified as: ", RecDevs_spec)
  } else collect_message("Recruitment Deviations is estimated for all dimensions")

  # Steepness
  if(input_list$data$rec_model == 0) {
    input_list$map$steepness_h <- factor(rep(NA, length(input_list$par$steepness_h)))
  } else if(!is.null(h_spec)) {
    if(rec_dd == 'global' && !h_spec %in% c("est_shared_r", "fix") && input_list$data$n_regions > 1) stop("Please specify a valid steepness option for global recruitment density dependence (should be est_shared_r or fix)!")
    # Share across regions and estimate
    if(h_spec == "est_shared_r") input_list$map$steepness_h <- factor(rep(1, length(input_list$par$steepness_h)))
    # Fix all steepness values
    if(h_spec == "fix") input_list$map$steepness_h <- factor(rep(NA, length(input_list$par$steepness_h)))
    if(!h_spec %in% c("est_shared_r", "fix"))  stop("Please specify a valid steepness option. These include: fix, est_shared_r. Conversely, leave at NULL to estimate all steepness values.")
    else collect_message("Steepness is specified as: ", h_spec)
  } else {
    if(input_list$data$rec_model == 1) {
      if(rec_dd == 'global' && input_list$data$n_regions > 1) stop("Please specify a valid steepness option for global recruitment density dependence (should be est_shared_r)!")
      input_list$map$steepness_h <- factor(c(1:input_list$data$n_regions))
    }
    collect_message("Steepness is estimated for all dimensions")
  }

  input_list$data$map_h_Pars <- as.numeric(input_list$map$steepness_h) # specify which ones are mapped off so they are only penalized once in priors

  # R0 proportions
  if(input_list$data$n_regions == 1) input_list$map$Rec_prop <- factor(rep(NA, length(input_list$par$Rec_prop)))

  # Print all messages if verbose is TRUE
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)

}

