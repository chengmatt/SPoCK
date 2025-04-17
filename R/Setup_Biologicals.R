#' Set up simulation stuff for biologicals
#'
#' @param n_yrs Number of years
#' @param n_regions Number of regions
#' @param n_sexes Number of sexes
#' @param n_sims Number of simulations
#' @param n_ages Number of ages
#' @param base_M_value Base natural mortality value dimensioned by regions, ages, sexes
#' @param M_pattern Natural mortality pattern. Options include: constant
#' @param base_WAA_values Base weight-at-age values dimensioned by regions, ages, sexes
#' @param WAA_pattern Weight-at-age pattern. Options include: constant
#' @param base_Maturity_AA_values Base maturity values dimensioned by regions, ages, sexes
#' @param Maturity_AA_pattern Maturity pattern. Options include: constant
#'
#' @export Setup_Sim_Biologicals
Setup_Sim_Biologicals <- function(n_sims = n_sims,
                                  n_yrs = n_yrs,
                                  n_regions = n_regions,
                                  n_ages = n_ages,
                                  n_sexes = n_sexes,
                                  base_M_value,
                                  M_pattern,
                                  base_WAA_values,
                                  WAA_pattern,
                                  base_Maturity_AA_values,
                                  Maturity_AA_pattern
                                  ) {

  # Create containers to store values
  M <- array(0, dim = c(n_yrs, n_regions, n_ages, n_sexes, n_sims))
  WAA <- array(0, dim = c(n_yrs, n_regions, n_ages, n_sexes, n_sims))
  Maturity_AA <- array(0, dim = c(n_yrs, n_regions, n_ages, n_sexes, n_sims))

  for(sim in 1:n_sims) {
    for(r in 1:n_regions) {
      for(y in 1:n_yrs) {
        for(s in 1:n_sexes) {

          # Natural mortality constant
          if(M_pattern == "constant") M[y,r,,s,sim] <- base_M_value[r,,s]

          # WAA constant
          if(WAA_pattern == "constant") WAA[y,r,,s,sim] <- base_WAA_values[r,,s]

          # Maturity constant
          if(Maturity_AA_pattern == "constant") Maturity_AA[y,r,,s,sim] <- base_Maturity_AA_values[r,,s]

        } # end s loop
      } # end y loop
    } # end r loop
  } # end sim loop

  # output into global environment
  M <<- M
  WAA <<- WAA
  Maturity_AA <<- Maturity_AA

}

#' Setup biological inputs for estimation model
#'
#' @param input_list List containing a data list, parameter list, and map list
#' @param WAA Inputs of weight at age dimensioned by n_region, n_year, n_ages, n_sexes
#' @param MatAA Inputs of maturity at age dimensioned by n_region, n_year, n_ages, n_sexes
#' @param AgeingError Inputs of ageing error matrix dimensioned by n_ages, n_ages. Default behavior uses an identity matrix (no ageing error)
#' @param Use_M_prior Whether or not to use natural mortality prior, == 0 don't use, == 1 use
#' @param M_prior Vector of natural mortality priors with the first element representing the mean in normal space, and the second element representing the sd.
#' @param fit_lengths Whether or not to fit length data, == 0 dont fit, == 1 fit
#' @param SizeAgeTrans Size age transition matrix dimensioned by n_regions, n_years, n_lens, n_ages, n_sexes
#' @param ... Parameter starting values for ln_M and M_offset if not using default
#' @param M_spec Character specifying options for how to estimate natural mortality. Default is NULL such that it is estimated for 2 sexes or only estimated for a single sex if it is a single sex model. Other options include "est_shared_s" which estimates the same natural mortality rate if n_sexes == 2. The other option is "fix" which fixes all natural mortality parameters.
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
                                  ...) {

  # Checking dimensions
  check_data_dimensions(WAA, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, what = 'WAA')
  check_data_dimensions(MatAA, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, what = 'MatAA')
  if(!is.null(AgeingError)) check_data_dimensions(AgeingError, n_ages = length(input_list$data$ages), what = 'AgeingError')
  if(fit_lengths == 1) check_data_dimensions(SizeAgeTrans, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_lens = length(input_list$data$lens), n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, what = 'SizeAgeTrans')

  input_list$data$WAA <- WAA
  input_list$data$MatAA <- MatAA
  if(is.null(AgeingError)) AgeingError <- diag(1, length(input_list$data$ages)) # if no inputs for ageing error, then create identity matrix
  input_list$data$AgeingError <- AgeingError
  input_list$data$fit_lengths <- fit_lengths
  input_list$data$SizeAgeTrans <- SizeAgeTrans
  input_list$data$Use_M_prior <- Use_M_prior
  input_list$data$M_prior <- M_prior

  if(!fit_lengths %in% c(0,1)) stop("Values for fit_lengths are not valid. They are == 0 (not used), or == 1 (used)")
  message("Length Composition data are: ", ifelse(fit_lengths == 0, "Not Used", "Used"))

  if(!Use_M_prior %in% c(0,1)) stop("Values for Use_M_prior are not valid. They are == 0 (don't use prior), or == 1 (use prior)")
  message("Natural Mortality priors are: ", ifelse(Use_M_prior == 0, "Not Used", "Used"))

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
    else message("Natural Mortality specified as: ", M_spec)
  } else if(input_list$data$n_sexes == 2) message("Natural Mortality is estiamted for both sexes") else message("Natural Mortality is estiamted")



  return(input_list)
}

