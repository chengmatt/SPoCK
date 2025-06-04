#' Set up simulation stuff for biologicals
#'
#' @param base_M_value Base natural mortality value dimensioned by regions, ages, sexes
#' @param M_pattern Natural mortality pattern. Options include: constant
#' @param base_WAA_values Base weight-at-age values dimensioned by regions, ages, sexes
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
                                  WAA_pattern,
                                  base_Maturity_AA_values,
                                  Maturity_AA_pattern,
                                  sim_list
                                  ) {

  # Create containers to store values
  M <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims))
  WAA <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims))
  Maturity_AA <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims))

  for(sim in 1:sim_list$n_sims) {
    for(r in 1:sim_list$n_regions) {
      for(y in 1:sim_list$n_yrs) {
        for(s in 1:sim_list$n_sexes) {

          # Natural mortality constant
          if(M_pattern == "constant") M[r,y,,s,sim] <- base_M_value[r,,s]

          # WAA constant
          if(WAA_pattern == "constant") WAA[r,y,,s,sim] <- base_WAA_values[r,,s]

          # Maturity constant
          if(Maturity_AA_pattern == "constant") Maturity_AA[r,y,,s,sim] <- base_Maturity_AA_values[r,,s]

        } # end s loop
      } # end y loop
    } # end r loop
  } # end sim loop

  # output into list
  sim_list$M <- M
  sim_list$WAA <- WAA
  sim_list$Maturity_AA <- Maturity_AA

  return(sim_list)

}

#' Setup biological inputs for estimation model
#'
#' @param input_list List containing a data list, parameter list, and map list
#' @param WAA Inputs of weight at age dimensioned by n_region, n_year, n_ages, n_sexes
#' @param MatAA Inputs of maturity at age dimensioned by n_region, n_year, n_ages, n_sexes
#' @param AgeingError Inputs of ageing error matrix dimensioned by number of modelled ages, number of observed composition ages Default behavior uses an identity matrix (no ageing error)
#' @param Use_M_prior Whether or not to use natural mortality prior, == 0 don't use, == 1 use
#' @param M_prior Vector of natural mortality priors with the first element representing the mean in normal space, and the second element representing the sd.
#' @param fit_lengths Whether or not to fit length data, == 0 dont fit, == 1 fit
#' @param SizeAgeTrans Size age transition matrix dimensioned by n_regions, n_years, n_lens, n_ages, n_sexes
#' @param M_spec Character specifying options for how to estimate natural mortality. Default is NULL such that it is estimated for 2 sexes or only estimated for a single sex if it is a single sex model. Other options include "est_shared_s" which estimates the same natural mortality rate if n_sexes == 2. The other option is "fix" which fixes all natural mortality parameters.
#' @param Fixed_natmort Fixed natural mortality array, dimensionsed by n_regions, n_yrs, n_ages, and n_sexes
#' @param ... Additional arguments for starting values of ln_M and M_offset (note that when M_spec is specified at 'fix', these will not be used. Instead, a user must supply a fixed natural mortality array using Fixed_natmort)
#'
#' @export Setup_Mod_Biologicals
Setup_Mod_Biologicals <- function(input_list,
                                  WAA,
                                  MatAA,
                                  AgeingError = NULL,
                                  Use_M_prior = 0,
                                  M_prior = NA,
                                  fit_lengths = 0,
                                  SizeAgeTrans = NA,
                                  M_spec = NULL,
                                  Fixed_natmort = NULL,
                                  ...
                                  ) {

  messages_list <<- character(0) # string to attach to for printing messages

  # Checking dimensions
  check_data_dimensions(WAA, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, what = 'WAA')
  check_data_dimensions(MatAA, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, what = 'MatAA')
  if(!is.null(AgeingError)) check_data_dimensions(AgeingError, n_ages = length(input_list$data$ages), what = 'AgeingError')
  if(fit_lengths == 1) check_data_dimensions(SizeAgeTrans, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_lens = length(input_list$data$lens), n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, what = 'SizeAgeTrans')
  if(!is.null(M_spec)) if(M_spec == 'fix') if(is.null(Fixed_natmort)) stop("Please provide a fixed natural mortality array dimensioned by n_regions, n_years, n_ages, and n_sexes!")
  if(!is.null(M_spec)) if(M_spec == 'fix') check_data_dimensions(Fixed_natmort, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, what = 'Fixed_natmort')

  input_list$data$WAA <- WAA
  input_list$data$MatAA <- MatAA
  if(is.null(AgeingError)) AgeingError <- diag(1, length(input_list$data$ages)) # if no inputs for ageing error, then create identity matrix
  input_list$data$AgeingError <- AgeingError
  input_list$data$fit_lengths <- fit_lengths
  input_list$data$SizeAgeTrans <- SizeAgeTrans
  input_list$data$Use_M_prior <- Use_M_prior
  input_list$data$M_prior <- M_prior
  input_list$data$Fixed_natmort <- Fixed_natmort

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

