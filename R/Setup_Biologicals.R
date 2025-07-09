#' Set up simulation stuff for biologicals
#'
#' @param base_M_value Base natural mortality value dimensioned by regions, ages, sexes
#' @param M_pattern Natural mortality pattern. Options include: constant
#' @param base_WAA_values Base weight-at-age values dimensioned by regions, ages, sexes
#' @param base_WAA_fish_values Base weight-at-age values for the fishery dimensioned by regions, ages, sexes, fishery fleets
#' @param WAA_pattern Weight-at-age pattern. Options include: constant
#' @param base_Maturity_AA_values Base maturity values dimensioned by regions, ages, sexes
#' @param Maturity_AA_pattern Maturity pattern. Options include: constant
#' @param sim_list Simulation list objects
#'
#' @export Setup_Sim_Biologicals
Setup_Sim_Biologicals <- function(
                                  base_M_value,
                                  M_pattern,
                                  base_WAA_values,
                                  base_WAA_fish_values,
                                  WAA_pattern,
                                  base_Maturity_AA_values,
                                  Maturity_AA_pattern,
                                  sim_list
                                  ) {

  # Create containers to store values
  M <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims))
  WAA <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims))
  WAA_fish <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_ages, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims))
  Maturity_AA <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims))

  for(sim in 1:sim_list$n_sims) {
    for(r in 1:sim_list$n_regions) {
      for(y in 1:sim_list$n_yrs) {
        for(s in 1:sim_list$n_sexes) {

          # Natural mortality constant
          if(M_pattern == "constant") M[r,y,,s,sim] <- base_M_value[r,,s]

          # WAA constant
          if(WAA_pattern == "constant") WAA[r,y,,s,sim] <- base_WAA_values[r,,s]
          if(WAA_pattern == "constant") for(f in 1:sim_list$n_fish_fleets) WAA_fish[r,y,,s,f,sim] <- base_WAA_fish_values[r,,s,f]

          # Maturity constant
          if(Maturity_AA_pattern == "constant") Maturity_AA[r,y,,s,sim] <- base_Maturity_AA_values[r,,s]

        } # end s loop
      } # end y loop
    } # end r loop
  } # end sim loop

  # output into list
  sim_list$M <- M
  sim_list$WAA <- WAA
  sim_list$WAA_fish <- WAA_fish
  sim_list$Maturity_AA <- Maturity_AA

  return(sim_list)

}

#' Setup biological inputs for estimation model
#'
#' @param input_list List containing data, parameter, and map lists for the model.
#' @param WAA Numeric array of weight-at-age (spawning), dimensioned \code{[n_regions, n_years, n_ages, n_sexes]}.
#' @param MatAA Numeric array of maturity-at-age, dimensioned \code{[n_regions, n_years, n_ages, n_sexes]}.
#' @param AgeingError Numeric matrix representing the ageing error transition matrix, dimensioned by \code{[number of modeled ages, number of observed composition ages]}. Defaults to identity matrix if not specified (no ageing error).
#' @param Use_M_prior Integer flag indicating whether to apply a natural mortality prior (\code{0} = no, \code{1} = yes).
#' @param M_prior Numeric vector of length two giving the mean (in normal space) and standard deviation of the natural mortality prior.
#' @param fit_lengths Integer flag indicating whether to fit length data (\code{0} = no, \code{1} = yes).
#' @param SizeAgeTrans Numeric array of size-at-age transition probabilities, dimensioned \code{[n_regions, n_years, n_lens, n_ages, n_sexes]}.
#' @param M_spec Character string specifying natural mortality estimation approach. Defaults to \code{NULL}, which estimates mortality for each sex independently. Other options:
#' \itemize{
#'   \item \code{"est_shared_s"}: Estimate a single natural mortality rate shared across sexes (if \code{n_sexes == 2}).
#'   \item \code{"fix"}: Fix all natural mortality parameters using the provided array.
#' }
#' @param Fixed_natmort Numeric array of fixed natural mortality values, dimensioned \code{[n_regions, n_years, n_ages, n_sexes]}. Required if \code{M_spec = "fix"}.
#' @param ... Additional arguments for starting values such as \code{ln_M} and \code{M_offset}. These are ignored if \code{M_spec = "fix"}.
#' @param Selex_Type Character string specifying whether selectivity is age or length-based. Default is age-based
#' \itemize{
#'   \item \code{"length"}: Length-based selectivity.
#'   \item \code{"age"}: Age-based selectivity
#' }
#' @param WAA_fish Numeric array of weight-at-age (fishery), dimensioned \code{[n_regions, n_years, n_ages, n_sexes, n_fish_fleets]}.
#' @param WAA_srv Numeric array of weight-at-age (survey), dimensioned \code{[n_regions, n_years, n_ages, n_sexes, n_srv_fleets]}.
#' @param addtocomp Numeric value for a constant to add to composition data. Default is 1e-3.
#'
#' @export Setup_Mod_Biologicals
Setup_Mod_Biologicals <- function(input_list,
                                  WAA,
                                  WAA_fish = NULL,
                                  WAA_srv = NULL,
                                  MatAA,
                                  addtocomp = 1e-3,
                                  AgeingError = NULL,
                                  Use_M_prior = 0,
                                  M_prior = NA,
                                  fit_lengths = 0,
                                  SizeAgeTrans = NA,
                                  Selex_Type = 'age',
                                  M_spec = NULL,
                                  Fixed_natmort = NULL,
                                  ...
                                  ) {

  messages_list <<- character(0) # string to attach to for printing messages

  # Checking dimensions
  check_data_dimensions(WAA, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, what = 'WAA')
  if(!is.null(WAA_fish)) check_data_dimensions(WAA_fish, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'WAA_fish')
  if(!is.null(WAA_srv)) check_data_dimensions(WAA_srv, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, n_srv_fleets = input_list$data$n_srv_fleets, what = 'WAA_srv')
  check_data_dimensions(MatAA, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, what = 'MatAA')
  if(!is.null(AgeingError)) check_data_dimensions(AgeingError, n_ages = length(input_list$data$ages), what = 'AgeingError')
  if(fit_lengths == 1) check_data_dimensions(SizeAgeTrans, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_lens = length(input_list$data$lens), n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, what = 'SizeAgeTrans')
  if(!is.null(M_spec)) if(M_spec == 'fix') if(is.null(Fixed_natmort)) stop("Please provide a fixed natural mortality array dimensioned by n_regions, n_years, n_ages, and n_sexes!")
  if(!is.null(M_spec)) if(M_spec == 'fix') check_data_dimensions(Fixed_natmort, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, what = 'Fixed_natmort')

  # Whether selectivity is age or length-based
  if(Selex_Type == 'age') {
    Selex_Type <- 0
    collect_message("Selectivity is aged-based.")
  } # if age based

  if(Selex_Type == 'length') {
    if(fit_lengths == 0) stop("Length composition data are not fit, but selectivity is length-based. This is not allowed. Please change to a valid option (either fit lengths or use age-based selectivity).")
    Selex_Type <- 1
    collect_message("Selectivity is length-based")
  } # if length based

  # setup fishery and survey specific weight at age (if not specified - just uses the WAA (spawning) already supplied)
  if(is.null(WAA_fish)) { # if no fishery WAA provided, use spawning WAA supplied
    WAA_fish <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), length(input_list$data$ages), input_list$data$n_sexes, input_list$data$n_fish_fleets))
    for(f in 1:input_list$data$n_fish_fleets) WAA_fish[,,,,f] <- WAA
    collect_message("WAA_fish was specified at NULL. Using the spawning WAA for WAA_fish")
  }

  # if no survey WAA provided, use spawning WAA supplied
  if(is.null(WAA_srv)) {
    WAA_srv <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), length(input_list$data$ages), input_list$data$n_sexes, input_list$data$n_srv_fleets))
    for(f in 1:input_list$data$n_srv_fleets) WAA_srv[,,,,f] <- WAA
    collect_message("WAA_srv was specified at NULL. Using the spawning WAA for WAA_srv")
  }

  input_list$data$WAA <- WAA
  input_list$data$WAA_fish <- WAA_fish
  input_list$data$WAA_srv <- WAA_srv
  input_list$data$MatAA <- MatAA
  if(is.null(AgeingError)) AgeingError <- diag(1, length(input_list$data$ages)) # if no inputs for ageing error, then create identity matrix
  input_list$data$AgeingError <- AgeingError
  input_list$data$fit_lengths <- fit_lengths
  input_list$data$SizeAgeTrans <- SizeAgeTrans
  input_list$data$Use_M_prior <- Use_M_prior
  input_list$data$M_prior <- M_prior
  input_list$data$Fixed_natmort <- Fixed_natmort
  input_list$data$Selex_Type <- Selex_Type
  input_list$data$addtocomp <- addtocomp

  # Input indicator for estimating or not estimating M
  if(is.null(M_spec) || M_spec == "est_ln_M_only") input_list$data$use_fixed_natmort <- 0
  else if(M_spec == "fix") input_list$data$use_fixed_natmort <- 1

  if(!fit_lengths %in% c(0,1)) stop("Values for fit_lengths are not valid. They are == 0 (not used), or == 1 (used)")
  collect_message("Length Composition data are: ", ifelse(fit_lengths == 0, "Not Used", "Used"))

  if(!Use_M_prior %in% c(0,1)) stop("Values for Use_M_prior are not valid. They are == 0 (don't use prior), or == 1 (use prior)")
  collect_message("Natural Mortality priors are: ", ifelse(Use_M_prior == 0, "Not Used", "Used"))

  if(fit_lengths == 1 & is.na(sum(SizeAgeTrans))) stop("Length composition are fit to, but the size-age transition matrix is NA")

  # Set up parameter input list
  starting_values <- list(...)
  if("ln_M" %in% names(starting_values)) input_list$par$ln_M <- starting_values$ln_M
  else input_list$par$ln_M <- log(0.5)

  if("M_offset" %in% names(starting_values)) input_list$par$M_offset <- starting_values$M_offset
  else input_list$par$M_offset <- 0

  # Setup mapping list
  if(input_list$data$n_sexes == 1) input_list$map$M_offset <- factor(NA) # fix sex-specific natural mortality offset if single sex
  if(!is.null(M_spec)) {
    # Estimate only 1 single natural mortality
    if(M_spec == 'est_ln_M_only') input_list$map$M_offset <- factor(NA)
    # Fix natural mortality rates for base value (ln_M) and offset (M_offset)
    if(M_spec == "fix") {
      input_list$map$ln_M <- factor(NA)
      input_list$map$M_offset <- factor(NA)
    }

    if(!M_spec %in% c('est_ln_M_only', 'fix')) stop("M_spec needs to be specified as either est_ln_M_only (only for a single sex) or fix")
    else collect_message("Natural Mortality specified as: ", M_spec)
  } else if(input_list$data$n_sexes == 2) collect_message("Natural Mortality is estiamted for both sexes") else collect_message("Natural Mortality is estiamted")

  # Print all messages if verbose is TRUE
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}

