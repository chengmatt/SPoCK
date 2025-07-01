#' Setup values and dimensions of fishing mortality
#'
#' @param sigmaC Observation error for catch
#' @param Fmort_pattern Fishing mortality pattern as a matrix dimensioned by region and fleet (constant, linear, one-way, two-way)
#' @param Fmort_start Fishing mortality start values as a matrix dimensioned by region and fleet
#' @param Fmort_fct Fishing mortality values for factor increases or decreases as a matrix dimensioned by region and fleet
#' @param proc_error If we want to add process error to fishing mortality
#' @param proc_error_sd value of logrnomal sd for process error to fishing mortality
#' @param init_F_vals Initial F values dimensioend by region and fleet
#' @param sim_list Simulation list objects
#'
#' @export Setup_Sim_FishMort
#' @importFrom stats rnorm
Setup_Sim_FishMort <- function(sim_list,
                               sigmaC,
                               init_F_vals,
                               Fmort_pattern,
                               Fmort_start,
                               Fmort_fct,
                               proc_error,
                               proc_error_sd
                                ) {

  # Fishing mortality stuff
  init_F <- array(0, dim = c(sim_list$n_regions, 1, sim_list$n_fish_fleets, sim_list$n_sims))
  Fmort <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs + 1, sim_list$n_fish_fleets, sim_list$n_sims))

  if(sim_list$run_feedback == TRUE) f_yrs <- sim_list$feedback_start_yr
  else f_yrs <- sim_list$n_yrs

  for(r in 1:sim_list$n_regions) {
    for(f in 1:sim_list$n_fish_fleets) {
      for(sim in 1:sim_list$n_sims) {

        # Set up initial F values
        init_F[r,1,f,sim] <- init_F_vals[r,f]

        if(Fmort_pattern[r,f] == 'constant') {
          Fmort[r,1:f_yrs,f,sim] <- Fmort_start[r,f]
        } # end if Fmort is constant

        if(Fmort_pattern[r,f] == "linear") {
          Fmort[r,1:f_yrs,f,sim] <- seq(Fmort_start[r,f], Fmort_start[r,f] * Fmort_fct[r,f], length.out = f_yrs)
        } # end if Fmort is linear increase of decrease

        if(Fmort_pattern[r,f] == "one-way") {
          Fmort[r,1:f_yrs,f,sim] <- c(seq(Fmort_start[r,f], Fmort_start[r,f] * Fmort_fct[r,f], length.out = round(f_yrs / 2)),
                               rep(Fmort_start[r,f] * Fmort_fct[r,f], ceiling(f_yrs / 2)))
        } # end if Fmort is a one way trip

        if(Fmort_pattern[r,f] == "two-way") {
          Fmort[r,1:f_yrs,f,sim] <- c(seq(Fmort_start[r,f], Fmort_start[r,f] * Fmort_fct[r,f], length.out = round(f_yrs / 2)),
                               seq(Fmort_start[r,f] * Fmort_fct[r,f], max(Fmort_start[r,f], Fmort_start[r,f] * 0.75), length.out = ceiling(f_yrs / 2)))
        } # end if Fmort is a two way trip

      } # end sim loop
    } # end f loop
  } # end r loop

  if(proc_error == TRUE) Fmort <- Fmort * exp(stats::rnorm(prod(dim(Fmort)), 0, proc_error_sd))

  # output variables into list
  sim_list$sigmaC <- sigmaC # Observation sd for catch
  sim_list$Fmort <- Fmort # fishing mortlaity pattern
  sim_list$init_F <- init_F # initial F value

  return(sim_list)

}

#' Setup fishery selectivity
#'
#' @param sel_model Fishery selectivity model dimensioned by region and fleet. Options include: logistic
#' @param fixed_fish_sel_pars Fixed parameters of fishery selectivity, dimensioned by region, sex, fishery fleet, and the max number of parameters needed for
#' for a defined fishery selectivity functional form out of all defined functional forms for the fishery
#' @param sim_list Simulation list
#'
#' @export Setup_Sim_FishSel
#'
Setup_Sim_FishSel <- function(sim_list,
                              sel_model,
                              fixed_fish_sel_pars
                              ) {

  # create fishery selectivity container
  fish_sel <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_ages, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims))

  for(sim in 1:sim_list$n_sims) {
    for(r in 1:sim_list$n_regions) {
      for(y in 1:sim_list$n_yrs) {
        for(f in 1:sim_list$n_fish_fleets) {
          for(s in 1:sim_list$n_sexes) {

            if(sel_model[r,f] == 'logistic') {
              a50 <- fixed_fish_sel_pars[r,s,f,1] # get a50
              k <- fixed_fish_sel_pars[r,s,f,2] # get k
              fish_sel[r,y,,s,f,sim] <- 1 / (1 + exp(-k * (1:sim_list$n_ages - a50)))
            } # end if logistic

            if(sel_model[r,f] == 'gamma') {
              amax <- fixed_fish_sel_pars[r,s,f,1] # get amax
              delta <- fixed_fish_sel_pars[r,s,f,2] # get delta
              p <- 0.5 * (sqrt( amax^2 + (4 * delta^2)) - amax)
              fish_sel[r,y,,s,f,sim] <- (1:n_ages / amax)^(amax/p) * exp( (amax - 1:sim_list$n_ages) / p )
            } # end if gamma

          }
        } # end f loop
      } # end y loop
    } # end r loop
  } # end sim loop


  # output fishery selectiviy
  sim_list$fish_sel <- fish_sel

  return(sim_list)

}


#' Setup fishing mortality and catch observations
#'
#' @param input_list A list containing data, parameters, and map lists used by the model.
#'
#' @param ObsCatch Numeric array of observed catches, dimensioned \code{[n_regions, n_years, n_fish_fleets]}.
#'
#' @param Catch_Type Integer matrix with dimensions \code{[n_years, n_fish_fleets]}, specifying catch data types:
#' \itemize{
#'   \item \code{0}: Use aggregated catch data for the year.
#'   \item \code{1}: Use region-specific catch data for the year.
#' }
#'
#' @param UseCatch Indicator array \code{[n_regions, n_years, n_fish_fleets]} specifying whether to include catch data in the fit:
#' \itemize{
#'   \item \code{0}: Do not use catch data.
#'   \item \code{1}: Use catch data and fit.
#' }
#'
#' @param Use_F_pen Integer flag indicating whether to apply a fishing mortality penalty:
#' \itemize{
#'   \item \code{0}: Do not apply penalty.
#'   \item \code{1}: Apply penalty.
#' }
#'
#' @param est_all_regional_F Integer flag indicating whether all regional fishing mortality deviations are estimated:
#' \itemize{
#'   \item \code{0}: Some fishing mortality deviations are aggregated across regions.
#'   \item \code{1}: All fishing mortality deviations are regional.
#' }
#'
#' @param Catch_Constant Numeric vector of length \code{n_fish_fleets} specifying constants to add to catch observations.
#'
#' @param sigmaC_spec Character string specifying observation error structure for catch data. Default behavior fixes \code{sigmaC} at a starting value of \code{1e-3} (log-scale \code{ln_sigmaC = log(1e-3)}) for all regions and fleets. Other options include:
#' \itemize{
#'   \item \code{"est_shared_f"}: Estimate \code{sigmaC} shared across fishery fleets.
#'   \item \code{"est_shared_r"}: Estimate \code{sigmaC} shared across regions but unique by fleet.
#'   \item \code{"est_shared_r_f"}: Estimate \code{sigmaC} shared across regions and fleets.
#'   \item \code{"fix"}: Fix \code{sigmaC} at the starting value.
#'   \item \code{"est_all"}: Estimate separate \code{sigmaC} for each region and fleet.
#' }
#'
#' @param sigmaF_spec Character string specifying process error structure for fishing mortality. Default fixes \code{sigmaF} at \code{1} on the log scale (i.e., \code{ln_sigmaF = 0}). Other options include:
#' \itemize{
#'   \item \code{"est_shared_f"}: Estimate \code{sigmaF} shared across fishery fleets.
#'   \item \code{"est_shared_r"}: Estimate \code{sigmaF} shared across regions but unique by fleet.
#'   \item \code{"est_shared_r_f"}: Estimate \code{sigmaF} shared across regions and fleets.
#'   \item \code{"fix"}: Fix \code{sigmaF} at the starting value.
#'   \item \code{"est_all"}: Estimate separate \code{sigmaF} for each region and fleet.
#' }
#'
#' @param sigmaF_agg_spec Character string specifying process error structure for aggregated fishing mortality. Default fixes \code{sigmaF_agg} at the starting value (log-scale \code{ln_sigmaF_agg}). Other options include:
#' \itemize{
#'   \item \code{"est_shared_f"}: Estimate \code{sigmaF_agg} shared across fishery fleets.
#'   \item \code{"fix"}: Fix at the starting value.
#'   \item \code{"est_all"}: Estimate separate parameters for each fishery fleet.
#' }
#'
#' @param ... Additional arguments specifying starting values for \code{ln_sigmaC}, \code{ln_sigmaF}, and \code{ln_sigmaF_agg}.
#'
#' @export Setup_Mod_Catch_and_F
#'
Setup_Mod_Catch_and_F <- function(input_list,
                                  ObsCatch,
                                  Catch_Type,
                                  UseCatch,
                                  Use_F_pen = 1,
                                  est_all_regional_F = 1,
                                  Catch_Constant = NULL,
                                  sigmaC_spec = "fix",
                                  sigmaF_spec = "fix",
                                  sigmaF_agg_spec = "fix",
                                  ...
                                  ) {

  messages_list <<- character(0) # string to attach to for printing messages

  # Dimension checking
  check_data_dimensions(ObsCatch, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsCatch')
  check_data_dimensions(Catch_Type, n_years = length(input_list$data$years), n_fish_fleets = input_list$data$n_fish_fleets, what = 'Catch_Type')
  check_data_dimensions(UseCatch, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseCatch')
  if(!is.null(Catch_Constant)) check_data_dimensions(Catch_Constant, n_fish_fleets = input_list$data$n_fish_fleets, what = 'Catch_Constant')

  else {
    if(est_all_regional_F == 0 && any(unique(Catch_Type) == 0)) collect_message("Catch is aggregated by region in some years, with a separate aggregated ln_F_Mean and ln_F_devs estimated in those years")
    if(est_all_regional_F == 1 && any(unique(Catch_Type) == 0)) collect_message("Catch is aggregated by region in some years, with a region specific ln_F_Mean and ln_F_devs estiamted, where these fishing mortalities are estimated using information from data (age and indices) in subsequent years")
    if(est_all_regional_F == 1 && any(unique(Catch_Type) == 1)) collect_message("Catch is region specific, with region specific ln_F_Mean and ln_F_devs")
  }

  if(!est_all_regional_F %in% c(0,1)) stop("est_all_regional_F incorrectly specified. Either set at 0 (not all regional Fs are estiamted) or 1 (all regional Fs are estimated)")
  else collect_message("Fishing mortality is estimated: ", ifelse(est_all_regional_F == 0, 'Not For All Regions', "For All Regions"))

  if(!Use_F_pen %in% c(0,1)) stop("Use_F_pen incorrectly specified. Either set at 0 (don't use F penalty) or 1 (use F penalty)")
  else collect_message("Fishing mortality penalty is: ", ifelse(Use_F_pen == 0, 'Not Used', "Used"))

  if(any(UseCatch == 0)) collect_message("User specified catch for some years and fleets to not be fit to, and ln_F_devs will not be estimated for those dimensions")

  # Input data list
  input_list$data$ObsCatch <- ObsCatch
  input_list$data$Catch_Type <- Catch_Type
  input_list$data$UseCatch <- UseCatch
  input_list$data$est_all_regional_F <- est_all_regional_F
  if(is.null(Catch_Constant)) Catch_Constant <- rep(0, input_list$data$n_fish_fleets)
  input_list$data$Catch_Constant <- Catch_Constant
  input_list$data$Use_F_pen <- Use_F_pen

  # Input parameters
  starting_values <- list(...)

  # Catch observation error
  if("ln_sigmaC" %in% names(starting_values)) input_list$par$ln_sigmaC <- starting_values$ln_sigmaC
  else input_list$par$ln_sigmaC <- array(log(1e-3), dim = c(input_list$data$n_regions, input_list$data$n_fish_fleets))

  # Process error fishing deviations for regional catch
  if("ln_sigmaF" %in% names(starting_values)) input_list$par$ln_sigmaF <- starting_values$ln_sigmaF
  else input_list$par$ln_sigmaF <- array(log(1), dim = c(input_list$data$n_regions, input_list$data$n_fish_fleets))

  # Process error fishing deviations for regional catch
  if("ln_sigmaF_agg" %in% names(starting_values)) input_list$par$ln_sigmaF <- starting_values$ln_sigmaF
  else input_list$par$ln_sigmaF_agg <- rep(log(1), input_list$data$n_fish_fleets)

  # Log mean fishing mortality
  if("ln_F_mean" %in% names(starting_values)) input_list$par$ln_F_mean <- starting_values$ln_F_mean
  else input_list$par$ln_F_mean <- array(log(0.1), dim = c(input_list$data$n_regions, input_list$data$n_fish_fleets))

  # Log fishing deviations
  if("ln_F_devs" %in% names(starting_values)) input_list$par$ln_F_devs <- starting_values$ln_F_devs
  else input_list$par$ln_F_devs <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets))

  # Log mean fishing mortality for aggregated catch
  if("ln_F_mean_AggCatch" %in% names(starting_values)) input_list$par$ln_F_mean_AggCatch <- starting_values$ln_F_mean_AggCatch
  else input_list$par$ln_F_mean_AggCatch <- rep(0, dim(input_list$data$Catch_Type[rowSums(input_list$data$Catch_Type) == 0, , drop = FALSE])[2])

  # Log fishing deviations for aggregated catch
  if("ln_F_devs_AggCatch" %in% names(starting_values)) input_list$par$ln_F_devs_AggCatch <- starting_values$ln_F_devs_AggCatch
  else input_list$par$ln_F_devs_AggCatch <- array(0, dim = dim(input_list$data$Catch_Type[rowSums(input_list$data$Catch_Type) == 0, , drop = FALSE]))

  # Setup mapping list

  # Observation error for catch
  map_sigmaC <- input_list$par$ln_sigmaC # initialize
  # Same sigmaC across fleets, but unique across regions
  if(sigmaC_spec == "est_shared_f") {
    map_sigmaC[1:input_list$data$n_regions,] <- 1:input_list$data$n_regions
    input_list$map$ln_sigmaC <- factor(map_sigmaC)
  }
  # Same sigmaC across regions, but unique across fleets
  if(sigmaC_spec == "est_shared_r") {
    map_sigmaC[,1:input_list$data$n_fish_fleets] <- 1:input_list$data$n_fish_fleets
    input_list$map$ln_sigmaC <- factor(map_sigmaC)
  }
  # Same sigmaC across regions and fleets
  if(sigmaC_spec == "est_shared_r_f") input_list$map$ln_sigmaC <- factor(rep(1, length(input_list$par$ln_sigmaC)))
  # Fixing sigmaC
  if(sigmaC_spec == "fix") input_list$map$ln_sigmaC <- factor(rep(NA, length(input_list$par$ln_sigmaC)))
  # Estimating all sigmaC
  if(sigmaC_spec == "est_all") input_list$map$ln_sigmaC <- factor(1:length(input_list$par$ln_sigmaC))
  if(!sigmaC_spec %in% c("est_shared_f", "est_shared_r", "est_shared_r_f", "fix", "est_all")) collect_message("sigmaC is specified as: ", sigmaC_spec)

  # Process error for fishing mortality
  map_sigmaF <- input_list$par$ln_sigmaF # initialize
  # Same sigmaF across fleets, but unique across regions
  if(sigmaF_spec == "est_shared_f") {
    map_sigmaF[1:input_list$data$n_regions,] <- 1:input_list$data$n_regions
    input_list$map$ln_sigmaF <- factor(map_sigmaF)
  }
  # Same sigmaF across regions, but unique across fleets
  if(sigmaF_spec == "est_shared_r") {
    map_sigmaF[,1:input_list$data$n_fish_fleets] <- 1:input_list$data$n_fish_fleets
    input_list$map$ln_sigmaF <- factor(map_sigmaF)
  }
  # Same sigmaF across regions and fleets
  if(sigmaF_spec == "est_shared_r_f") input_list$map$ln_sigmaF <- factor(rep(1, length(input_list$par$ln_sigmaF)))
  # Fixing sigmaF
  if(sigmaF_spec == "fix") input_list$map$ln_sigmaF <- factor(rep(NA, length(input_list$par$ln_sigmaF)))
  # Estimating all sigmaF
  if(sigmaF_spec == "est_all") input_list$map$ln_sigmaF <- factor(1:length(input_list$par$ln_sigmaF))
  collect_message("sigmaF is specified as: ", sigmaF_spec)

  # Process error for aggregated fishing mortality / catch options
  if(input_list$data$est_all_regional_F == 0) {
    # Same sigmaF across fleets
    if(sigmaF_agg_spec == "est_shared_f") input_list$map$ln_sigmaF <- factor(rep(1, length(input_list$par$ln_sigmaF_agg)))
    # Fixing ln_sigmaF_agg
    if(sigmaF_agg_spec == "fix") input_list$map$ln_sigmaF_agg <- factor(rep(NA, length(input_list$par$ln_sigmaF_agg)))
    # Estimating all ln_sigmaF_agg
    if(sigmaF_agg_spec == "est_all") input_list$map$ln_sigmaF_agg <- factor(1:length(input_list$par$ln_sigmaF_agg))
  } # end if some fishing mortality deviations are aggregated

  # If all fishing mortality deviations are regional, then don't estimate this
  if(input_list$data$est_all_regional_F == 1 && sigmaF_agg_spec != 'fix') stop("Fishing mortality is specified to be aggregated for all regions, but sigmaF_agg_spec is not specified at `fix`!")
  if(input_list$data$est_all_regional_F == 1) input_list$map$ln_sigmaF_agg <- factor(rep(NA, length(input_list$par$ln_sigmaF_agg)))
  collect_message("sigmaF_agg is specified as: ", sigmaF_agg_spec)

  # Mapping for fishing mortality deviations
  F_dev_map <- input_list$par$ln_F_devs # initialize for mapping
  F_dev_map[] <- NA
  F_dev_counter <- 1
  for(r in 1:input_list$data$n_regions) {
    for(y in 1:length(input_list$data$years)) {
      for(f in 1:input_list$data$n_fish_fleets) {
        # if we are not using catch or its specified as aggregated catch, then turn off estimation of these f devs
        if(input_list$data$UseCatch[r,y,f] == 0) {
          F_dev_map[r,y,f] <- NA
        }
        # if we are using regional catch or its specified as regional catch, then estimate the f devs
        if(input_list$data$UseCatch[r,y,f] == 1) {
          F_dev_map[r,y,f] <- F_dev_counter
          F_dev_counter <- F_dev_counter + 1
        } # end if
      } # end f
    } # end y
  } # end r

  # If we have regional aggregated catch and some F devs are not regional
  if(any(input_list$dat$Catch_Type == 0) && input_list$data$est_all_regional_F == 0) {
    Catch_Type_map <- input_list$data$Catch_Type # initialize for mapping
    # Loop through to find years and fleets with aggregated catch
    for(f in 1:input_list$data$n_fish_fleets) {
      agg_tmp <- which(Catch_Type_map[,f] == 0) # years with aggregated catch for a given fleet
      F_dev_map[,agg_tmp,f] <- NA # specify as NA for 0s
      F_dev_map[,-agg_tmp,f] <- 1:length(F_dev_map[,-agg_tmp,f]) # specify as unique values for others
    } # end f loop
  }

  input_list$map$ln_F_devs <- factor(F_dev_map)

  # If we are estimating some aggregated F devs
  if(input_list$data$est_all_regional_F == 0) {
    input_list$map$ln_F_devs_AggCatch <- factor(1:length(input_list$par$ln_F_devs_AggCatch))
    input_list$map$ln_F_mean_AggCatch <- factor(1:input_list$data$n_fish_fleets)
  }

  # If we are estimating regional F devs for all
  if(input_list$data$est_all_regional_F == 1) {
    input_list$map$ln_F_devs_AggCatch <- factor(array(NA, dim = dim(input_list$data$Catch_Type[rowSums(input_list$data$Catch_Type) == 0, , drop = FALSE])))
    input_list$map$ln_F_mean_AggCatch <- factor(rep(NA, input_list$data$n_fish_fleets))
  }

  # Print all messages if verbose is TRUE
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}


#' Setup observed fishery indices and composition data (age and length comps)
#'
#' @param input_list List containing a data list, parameter list, and map list
#' @param ObsFishIdx Observed fishery index data as a numeric array with dimensions
#' \code{[n_regions, n_years, n_fish_fleets]}.
#'
#' @param ObsFishIdx_SE Standard errors associated with \code{ObsFishIdx},
#' also dimensioned \code{[n_regions, n_years, n_fish_fleets]}.
#'
#' @param UseFishIdx Logical or binary indicator array (\code{[n_regions, n_years, n_fish_fleets]})
#' specifying whether to include a fishery index in the likelihood (\code{1}) or ignore it (\code{0}).
#'
#' @param ObsFishAgeComps Observed fishery age composition data as a numeric array with dimensions
#' \code{[n_regions, n_years, n_ages, n_sexes, n_fish_fleets]}. Values should reflect counts or proportions
#' (not required to sum to 1, but should be on a comparable scale).
#'
#' @param UseFishAgeComps Indicator array (\code{[n_regions, n_years, n_fish_fleets]}) specifying whether
#' to fit fishery age composition data (\code{1}) or ignore it (\code{0}).
#'
#' @param ObsFishLenComps Observed fishery length composition data as a numeric array with dimensions
#' \code{[n_regions, n_years, n_lens, n_sexes, n_fish_fleets]}. Values should reflect counts or proportions.
#'
#' @param UseFishLenComps Indicator array (\code{[n_regions, n_years, n_fish_fleets]}) specifying whether
#' to fit fishery length composition data (\code{1}) or ignore it (\code{0}).
#'
#' @param FishAgeComps_LikeType Character vector of length \code{n_fish_fleets} specifying the likelihood
#' type used for fishery age composition data. Options include \code{"Multinomial"}, \code{"Dirichlet-Multinomial"},
#' and \code{"iid-Logistic-Normal"}. Use \code{"none"} to omit the likelihood.
#'
#' @param FishLenComps_LikeType Same as \code{FishAgeComps_LikeType}, but for fishery length composition data.
#'
#' @param FishAgeComps_Type Character vector specifying how age compositions are structured by fleet and year range.
#' Options include:
#' \itemize{
#'   \item \code{"agg"}: Aggregated across regions and sexes.
#'   \item \code{"spltRspltS"}: Split by region and by sex (compositions sum to 1 within region-sex group).
#'   \item \code{"spltRjntS"}: Split by region but summed jointly across sexes.
#'   \item \code{"none"}: No composition data used.
#' }
#' Format each element as \code{"<type>_Year_<start>-<end>_Fleet_<fleet number>"}
#' (e.g., \code{"agg_Year_1-10_Fleet_1"}).
#'
#' @param FishLenComps_Type Same as \code{FishAgeComps_Type}, but for length compositions.
#'
#' @param FishAge_comp_agg_type Optional integer vector of length \code{n_fish_fleets} specifying
#' the order of operations for aggregating age compositions when \code{FishAgeComps_Type == "agg"}.
#' \itemize{
#'   \item \code{0}: Normalize, then aggregate, then apply ageing error, then normalize again.
#'   \item \code{1}: Aggregate first, normalize, then apply ageing error.
#' }
#' Default is \code{NULL}.
#'
#' @param FishLen_comp_agg_type Optional integer vector of length \code{n_fish_fleets} specifying
#' the order of operations for aggregating length compositions.
#' \itemize{
#'   \item \code{0}: Do not normalize before applying size–age transition.
#'   \item \code{1}: Normalize before applying size–age transition.
#' }
#' Default is \code{NULL}.
#'
#' @param fish_idx_type Character vector of length \code{n_fish_fleets} specifying the type of index data.
#' Options are \code{"abd"} for abundance, \code{"biom"} for biomass, and \code{"none"} if no index is available.
#'
#' @param ISS_FishAgeComps Input sample size for age compositions, array dimensioned
#' \code{[n_regions, n_years, n_sexes, n_fish_fleets]}. Required if observed age comps are normalized
#' (i.e., sum to 1), to correctly scale the contribution to the likelihood.
#'
#' @param ISS_FishLenComps Same as \code{ISS_FishAgeComps}, but for length compositions.
#'
#' @param ... Additional arguments specifying starting values for overdispersion parameters
#' (e.g., \code{ln_FishAge_theta}, \code{ln_FishLen_theta}, \code{ln_FishAge_theta_agg}, \code{ln_FishLen_theta_agg}).
#'
#' @export Setup_Mod_FishIdx_and_Comps
#' @importFrom stringr str_detect
Setup_Mod_FishIdx_and_Comps <- function(input_list,
                                        ObsFishIdx,
                                        ObsFishIdx_SE,
                                        fish_idx_type,
                                        UseFishIdx,
                                        ObsFishAgeComps,
                                        UseFishAgeComps,
                                        ISS_FishAgeComps = NULL,
                                        ObsFishLenComps,
                                        UseFishLenComps,
                                        ISS_FishLenComps = NULL,
                                        FishAgeComps_LikeType,
                                        FishLenComps_LikeType,
                                        FishAgeComps_Type,
                                        FishLenComps_Type,
                                        FishAge_comp_agg_type = NULL,
                                        FishLen_comp_agg_type = NULL,
                                        ...
                                        ) {

  messages_list <<- character(0) # string to attach to for printing messages

  # Checking dimensions
  check_data_dimensions(ObsFishIdx, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsFishIdx')
  check_data_dimensions(ObsFishIdx_SE, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsFishIdx_SE')
  check_data_dimensions(UseFishIdx, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseFishIdx')
  check_data_dimensions(ObsFishAgeComps, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsFishAgeComps')
  check_data_dimensions(UseFishAgeComps, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseFishAgeComps')
  check_data_dimensions(UseFishLenComps, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseFishLenComps')
  check_data_dimensions(ObsFishLenComps, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_lens = length(input_list$data$lens), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsFishLenComps')
  if(!is.null(ISS_FishAgeComps)) check_data_dimensions(ISS_FishAgeComps, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ISS_FishAgeComps')
  if(!is.null(ISS_FishLenComps)) check_data_dimensions(ISS_FishLenComps, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ISS_FishLenComps')
  check_data_dimensions(fish_idx_type, n_fish_fleets = input_list$data$n_fish_fleets, what = 'fish_idx_type')
  check_data_dimensions(FishAgeComps_LikeType, n_fish_fleets = input_list$data$n_fish_fleets, what = 'FishAgeComps_LikeType')
  check_data_dimensions(FishLenComps_LikeType, n_fish_fleets = input_list$data$n_fish_fleets, what = 'FishLenComps_LikeType')

  # Checking specifications for index, age comps, length comps
  if(!all(fish_idx_type %in% c("biom", "abd", "none"))) stop("Invalid specification for fish_idx_type. Should be either abd, biom, or none")
  if(!all(FishAgeComps_LikeType %in% c("none", "Multinomial", "Dirichlet-Multinomial", "iid-Logistic-Normal", "1d-Logistic-Normal", "2d-Logistic-Normal", "3d-Logistic-Normal")))
    stop("Invalid specification for FishAgeComps_LikeType Should be either none, Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal, 1d-Logistic-Normal, 2d-Logistic-Normal, 3d-Logistic-Normal")
  if(!all(FishLenComps_LikeType %in% c("none", "Multinomial", "Dirichlet-Multinomial", "iid-Logistic-Normal", "1d-Logistic-Normal", "2d-Logistic-Normal", "3d-Logistic-Normal")))
    stop("Invalid specification for FishLenComps_LikeType Should be either none, Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal, 1d-Logistic-Normal, 2d-Logistic-Normal, 3d-Logistic-Normal")

  # Specifying fishery index type
  fish_idx_type_vals <- array(NA, dim = c(input_list$data$n_regions, input_list$data$n_fish_fleets))
  for(f in 1:ncol(fish_idx_type_vals)) {
    if(fish_idx_type[f] == 'biom') fish_idx_type_vals[,f] <- 1 # biomass
    if(fish_idx_type[f] == 'abd') fish_idx_type_vals[,f] <- 0 # abundance
    if(fish_idx_type[f] == 'none') fish_idx_type_vals[,f] <- 999 # none
    collect_message(paste("Fishery Index", "for fishery fleet", f, "specified as:" , fish_idx_type[f]))
  } # end f loop

  # Specifying fishery age composition values
  comp_fishage_like_vals <- vector()
  for(f in 1:input_list$data$n_fish_fleets) {
    if(FishAgeComps_LikeType[f] == 'none') comp_fishage_like_vals <- c(comp_fishage_like_vals, 999)
    if(FishAgeComps_LikeType[f] == "Multinomial") comp_fishage_like_vals <- c(comp_fishage_like_vals, 0)
    if(FishAgeComps_LikeType[f] == "Dirichlet-Multinomial") comp_fishage_like_vals <- c(comp_fishage_like_vals, 1)
    if(FishAgeComps_LikeType[f] == "iid-Logistic-Normal") comp_fishage_like_vals <- c(comp_fishage_like_vals, 2)
    if(FishAgeComps_LikeType[f] == "1d-Logistic-Normal") comp_fishage_like_vals <- c(comp_fishage_like_vals, 3)
    if(FishAgeComps_LikeType[f] == "2d-Logistic-Normal") comp_fishage_like_vals <- c(comp_fishage_like_vals, 4)
    if(FishAgeComps_LikeType[f] == "3d-Logistic-Normal") comp_fishage_like_vals <- c(comp_fishage_like_vals, 5)
    collect_message(paste("Fishery Age Composition Likelihoods", "for fishery fleet", f, "specified as:" , FishAgeComps_LikeType[f]))
  } # end f loop

  # Specifying fishery len composition values
  comp_fishlen_like_vals <- vector()
  for(f in 1:input_list$data$n_fish_fleets) {
    if(FishLenComps_LikeType[f] == 'none') comp_fishlen_like_vals <- c(comp_fishlen_like_vals, 999)
    if(FishLenComps_LikeType[f] == "Multinomial") comp_fishlen_like_vals <- c(comp_fishlen_like_vals, 0)
    if(FishLenComps_LikeType[f] == "Dirichlet-Multinomial") comp_fishlen_like_vals <- c(comp_fishlen_like_vals, 1)
    if(FishLenComps_LikeType[f] == "iid-Logistic-Normal") comp_fishlen_like_vals <- c(comp_fishlen_like_vals, 2)
    if(FishLenComps_LikeType[f] == "1d-Logistic-Normal") comp_fishlen_like_vals <- c(comp_fishlen_like_vals, 3)
    if(FishLenComps_LikeType[f] == "2d-Logistic-Normal") comp_fishlen_like_vals <- c(comp_fishlen_like_vals, 4)
    if(FishLenComps_LikeType[f] == "3d-Logistic-Normal") comp_fishlen_like_vals <- c(comp_fishlen_like_vals, 5)
    collect_message(paste("Fishery Length Composition Likelihoods", "for fishery fleet", f, "specified as:" , FishLenComps_LikeType[f]))
  } # end f loop

  # Specifying composition type (ages)
  FishAgeComps_Type_Mat <- array(NA, dim = c(length(input_list$data$years), input_list$data$n_fish_fleets))
  for(i in 1:length(FishAgeComps_Type)) {

    # Extract out components from list
    tmp <- FishAgeComps_Type[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    comps_type_tmp <- tmp_vec[1] # get composition type

    # get year ranges
    if(!str_detect(tmp, "terminal")) { # if not terminal year
      year_range <- as.numeric(unlist(strsplit(tmp_vec[3], "-")))
      years <- year_range[1]:year_range[2] # get sequence of years
    } else { # if terminal year
      year_range <- unlist(strsplit(tmp_vec[3], '-'))[1] # get year range
      years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
    }

    fleet <- as.numeric(tmp_vec[5]) # extract fleet index

    # Checking character string
    if(!comps_type_tmp %in% c("agg", "spltRspltS", "spltRjntS", 'none')) stop("FishAgeComps_Type not specified correctly. Must be one of: agg, spltRspltS, spltRjntS, none")
    if(!fleet %in% c(1:input_list$data$n_fish_fleets)) stop("Invalid fleet specified for FishAgeComps_Type. This needs to be specified as CompType_Year_x-y_Fleet_x")

    if(comps_type_tmp == "agg") comps_type_val <- 0
    if(comps_type_tmp == "spltRspltS") comps_type_val <- 1
    if(comps_type_tmp == "spltRjntS") comps_type_val <- 2
    if(comps_type_tmp == "none") comps_type_val <- 999

    # input into matrix
    FishAgeComps_Type_Mat[years,fleet] <- comps_type_val
  } # end i

  if(any(is.na(FishAgeComps_Type_Mat))) stop("FishAgeComps_Type is returning an NA. Did you update the year range of FishAgeComps_Type?")

  FishLenComps_Type_Mat <- array(NA, dim = c(length(input_list$data$years), input_list$data$n_fish_fleets))
  for(i in 1:length(FishLenComps_Type)) {

    # Extract out components from list
    tmp <- FishLenComps_Type[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    comps_type_tmp <- tmp_vec[1] # get composition type

    # get year ranges
    if(!str_detect(tmp, "terminal")) { # if not terminal year
      year_range <- as.numeric(unlist(strsplit(tmp_vec[3], "-")))
      years <- year_range[1]:year_range[2] # get sequence of years
    } else { # if terminal year
      year_range <- unlist(strsplit(tmp_vec[3], '-'))[1] # get year range
      years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
    }

    fleet <- as.numeric(tmp_vec[5]) # extract fleet index

    # define composition types
    if(!comps_type_tmp %in% c("agg", "spltRspltS", "spltRjntS", 'none')) stop("FishLenComps_Type not specified correctly. Must be one of: agg, spltRspltS, spltRjntS, none")
    if(!fleet %in% c(1:input_list$data$n_fish_fleets)) stop("Invalid fleet specified for FishLenComps_Type This needs to be specified as CompType_Year_x-y_Fleet_x")

    # define composition types
    if(comps_type_tmp == "agg") comps_type_val <- 0
    if(comps_type_tmp == "spltRspltS") comps_type_val <- 1
    if(comps_type_tmp == "spltRjntS") comps_type_val <- 2
    if(comps_type_tmp == "none") comps_type_val <- 999

    # input into matrix
    FishLenComps_Type_Mat[years,fleet] <- comps_type_val
  } # end i

  if(any(is.na(FishLenComps_Type_Mat))) stop("FishLenComps_Type_Mat is returning an NA. Did you update the year range of FishLenComps_Type_Mat?")

  # Input data list
  input_list$data$ObsFishIdx <- ObsFishIdx
  input_list$data$ObsFishIdx_SE <- ObsFishIdx_SE
  input_list$data$UseFishIdx <- UseFishIdx
  input_list$data$fish_idx_type <- fish_idx_type_vals
  input_list$data$ObsFishAgeComps <- ObsFishAgeComps
  input_list$data$UseFishAgeComps <- UseFishAgeComps
  input_list$data$ObsFishLenComps <- ObsFishLenComps
  input_list$data$UseFishLenComps <- UseFishLenComps
  input_list$data$FishAgeComps_LikeType <- comp_fishage_like_vals
  input_list$data$FishLenComps_LikeType <- comp_fishlen_like_vals
  input_list$data$FishAgeComps_Type <- FishAgeComps_Type_Mat
  input_list$data$FishLenComps_Type <- FishLenComps_Type_Mat

  # initialize how to aggregate fishery age comps
  if(is.null(FishAge_comp_agg_type)) {
    input_list$data$FishAge_comp_agg_type <- rep(NA, input_list$data$n_fish_fleets)
  } else input_list$data$FishAge_comp_agg_type <- FishAge_comp_agg_type

  # Initialize how to aggregate fishery length comps
  if(is.null(FishLen_comp_agg_type)) {
    input_list$data$FishLen_comp_agg_type <- rep(NA, input_list$data$n_fish_fleets)
  } else input_list$data$FishLen_comp_agg_type <- FishLen_comp_agg_type

  # Setup ISS stuff given differences in how composition data are structured
  if(is.null(ISS_FishAgeComps)) {
    collect_message("No ISS is specified for FishAgeComps. ISS weighting is calculated by summing up values from ObsFishAgeComps each year")
    ISS_FishAgeComps <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_sexes, input_list$data$n_fish_fleets))
    for(y in 1:length(input_list$data$years)) {
      for(f in 1:input_list$data$n_fish_fleets) {
        # Fishery Age Compositions
        # if aggregated across sexes and regions (0) or joint across sexes and regions (3)
        if(input_list$data$FishAgeComps_Type[y,f] %in% c(0, 3)) ISS_FishAgeComps[1,y,1,f] <- sum(input_list$data$ObsFishAgeComps[,y,,,f])
        # if split by region and sex
        if(input_list$data$FishAgeComps_Type[y,f] == 1) ISS_FishAgeComps[,y,,f] <- apply(input_list$data$ObsFishAgeComps[,y,,,f, drop = FALSE], c(1,4), sum)
        # if split by region, joint by sex
        if(input_list$data$FishAgeComps_Type[y,f] == 2) ISS_FishAgeComps[,y,1,f] <- apply(input_list$data$ObsFishAgeComps[,y,,,f, drop = FALSE], 1, sum)
      } # end f loop
    } # end y loop
  }

  if(is.null(ISS_FishLenComps)) {
    collect_message("No ISS is specified for FishLenComps. ISS weighting is calculated by summing up values from ObsFishLenComps each year")
    ISS_FishLenComps <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_sexes, input_list$data$n_fish_fleets))
    for(y in 1:length(input_list$data$years)) {
      for(f in 1:input_list$data$n_fish_fleets) {
        # Fishery Length Compositions
        # if aggregated across sexes and regions (0) or joint across sexes and regions (3)
        if(input_list$data$FishLenComps_Type[y,f] %in% c(0, 3)) ISS_FishLenComps[1,y,1,f] <- sum(input_list$data$ObsFishLenComps[,y,,,f])
        # if split by region and sex
        if(input_list$data$FishLenComps_Type[y,f] == 1) ISS_FishLenComps[,y,,f] <- apply(input_list$data$ObsFishLenComps[,y,,,f, drop = FALSE], c(1,4), sum)
        # if split by region, joint by sex
        if(input_list$data$FishLenComps_Type[y,f] == 2) ISS_FishLenComps[,y,1,f] <- apply(input_list$data$ObsFishLenComps[,y,,,f, drop = FALSE], 1, sum)
      } # end f loop
    } # end y loop
  }

  # Input into list
  input_list$data$ISS_FishAgeComps <- ISS_FishAgeComps
  input_list$data$ISS_FishLenComps <- ISS_FishLenComps

  # input parameter list
  starting_values <- list(...)

  # Dispersion parameters for the fishery age comps
  if("ln_FishAge_theta" %in% names(starting_values)) input_list$par$ln_FishAge_theta <- starting_values$ln_FishAge_theta
  else input_list$par$ln_FishAge_theta <- array(0, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets))

  # logistic normal correlation parameters for fishery age comps
  if("FishAge_corr_pars" %in% names(starting_values)) input_list$par$FishAge_corr_pars <- starting_values$FishAge_corr_pars
  else input_list$par$FishAge_corr_pars <- array(0.01, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets, 3))

  # aggregated
  if("ln_FishAge_theta_agg" %in% names(starting_values)) input_list$par$ln_FishAge_theta_agg <- starting_values$ln_FishAge_theta_agg
  else input_list$par$ln_FishAge_theta_agg <- array(0, dim = c(input_list$data$n_fish_fleets))

  # aggregated correlation parameters
  if("FishAge_corr_pars_agg" %in% names(starting_values)) input_list$par$FishAge_corr_pars_agg <- starting_values$FishAge_corr_pars_agg
  else input_list$par$FishAge_corr_pars_agg <- array(0.01, dim = c(input_list$data$n_fish_fleets))

  # Dispersion parameters for fishery length comps
  if("ln_FishLen_theta" %in% names(starting_values)) input_list$par$ln_FishLen_theta <- starting_values$ln_FishLen_theta
  else input_list$par$ln_FishLen_theta <- array(0, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets))

  # logistic normal correlation parameters for fishery length comps
  if("FishLen_corr_pars" %in% names(starting_values)) input_list$par$FishLen_corr_pars <- starting_values$FishLen_corr_pars
  else input_list$par$FishLen_corr_pars <- array(0.01, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets, 3))

  # aggregated
  if("ln_FishLen_theta_agg" %in% names(starting_values)) input_list$par$ln_FishLen_theta_agg <- starting_values$ln_FishLen_theta_agg
  else input_list$par$ln_FishLen_theta_agg <- array(0, dim = c(input_list$data$n_fish_fleets))

  if("FishLen_corr_pars_agg" %in% names(starting_values)) input_list$par$FishLen_corr_pars_agg <- starting_values$FishLen_corr_pars_agg
  else input_list$par$FishLen_corr_pars_agg <- array(0.01, dim = c(input_list$data$n_fish_fleets))

  # Setup counters
  counter_fishage_agg <- 1
  counter_fishage <- 1
  counter_fishlen_agg <- 1
  counter_fishlen <- 1
  counter_fishage_corr <- 1
  counter_fishage_corr_agg <- 1
  counter_fishlen_corr <- 1
  counter_fishlen_corr_agg <- 1

  # Setup mapping list
  map_FishAge_theta <- input_list$par$ln_FishAge_theta # initialize array to set up mapping
  map_FishLen_theta <- input_list$par$ln_FishLen_theta # initialize array to set up mapping
  map_FishAge_theta_agg <- input_list$par$ln_FishAge_theta_agg # initialize array to set up mapping
  map_FishLen_theta_agg <- input_list$par$ln_FishLen_theta_agg # initialize array to set up mapping
  map_FishAge_corr_pars <- input_list$par$FishAge_corr_pars # initialize array to set up mapping
  map_FishAge_corr_pars_agg <- input_list$par$FishAge_corr_pars_agg # initialize array to set up mapping
  map_FishLen_corr_pars <- input_list$par$FishLen_corr_pars # initialize array to set up mapping
  map_FishLen_corr_pars_agg <- input_list$par$FishLen_corr_pars_agg # initialize array to set up mapping

  # Initialize these as NAs
  map_FishAge_theta[] <- NA; map_FishLen_theta[] <- NA; map_FishAge_theta_agg[] <- NA; map_FishLen_theta_agg[] <- NA; map_FishAge_corr_pars[] <- NA; map_FishAge_corr_pars_agg[] <- NA; map_FishLen_corr_pars[] <- NA; map_FishLen_corr_pars_agg[] <- NA

  for(f in 1:input_list$data$n_fish_fleets) {
    # If we are using a multinomial or there aren't any age comps for a given fleet
    if(input_list$data$FishAgeComps_LikeType[f] == 0 || sum(input_list$data$UseFishAgeComps[,,f]) == 0) {
      map_FishAge_theta[,,f] <- NA
      map_FishAge_theta_agg[f] <- NA
      map_FishAge_corr_pars[,,f,] <- NA
      map_FishAge_corr_pars_agg[f] <- NA
    }

    # If we are using a multinomial or there aren't any length comps for a given fleet
    if(input_list$data$FishLenComps_LikeType[f] == 0 || sum(input_list$data$UseFishLenComps[,,f]) == 0) {
      map_FishLen_theta[,,f] <- NA
      map_FishLen_theta_agg[f] <- NA
      map_FishLen_corr_pars[,,f,] <- NA
      map_FishLen_corr_pars_agg[f] <- NA
    }

    # get unique fishery age comp types
    fishage_comp_type <- unique(input_list$data$FishAgeComps_Type[,f])
    fishlen_comp_type <- unique(input_list$data$FishLenComps_Type[,f])

    # If aggregated (ages)
    if(any(fishage_comp_type == 0) && input_list$data$FishAgeComps_LikeType[f] != 0) {
      map_FishAge_theta_agg[f] <- counter_fishage_agg
      counter_fishage_agg <- counter_fishage_agg + 1 # aggregated

      # correlation parameters if any aggregated fishery ages
      if(input_list$data$FishAgeComps_LikeType[f] == 3) {
        map_FishAge_corr_pars_agg[f] <- counter_fishage_corr_agg
        counter_fishage_corr_agg <- counter_fishage_corr_agg + 1 # aggregated
      }
    }

    # If aggregated (lengths)
    if(any(fishlen_comp_type == 0) && input_list$data$FishLenComps_LikeType[f] != 0) {
      map_FishLen_theta_agg[f] <- counter_fishlen_agg
      counter_fishlen_agg <- counter_fishlen_agg + 1 # aggregated

      # correlation parameters if any aggregated fishery lengths
      if(input_list$data$FishLenComps_LikeType[f] == 3) {
        map_FishLen_corr_pars_agg[f] <- counter_fishlen_corr_agg
        counter_fishlen_corr_agg <- counter_fishlen_corr_agg + 1 # aggregated
      }
    }

    # Loop through to make sure mapping stuff off correctly
    for(r in 1:input_list$data$n_regions) {
      for(s in 1:input_list$data$n_sexes) {

        # if split by sex and region (ages)
        if(any(fishage_comp_type == 1) && input_list$data$FishAgeComps_LikeType[f] != 0) {
          map_FishAge_theta[r,s,f] <- counter_fishage
          counter_fishage <- counter_fishage + 1 # split by sex and region

          # if logistic normal 1dar1 by age only (unique age correlations for age each region and sex)
          if(input_list$data$FishAgeComps_LikeType[f] == 3) {
            map_FishAge_corr_pars[r,s,f,1] <- counter_fishage_corr
            counter_fishage_corr <- counter_fishage_corr + 1
          }
        }

        # if split by sex and region (lengths)
        if(any(fishlen_comp_type == 1) && input_list$data$FishLenComps_LikeType[f] != 0) {
          map_FishLen_theta[r,s,f] <- counter_fishlen
          counter_fishlen <- counter_fishlen + 1 # split by sex and region

          # if logistic normal 1dar1 by length only (unique length correlations for length each region and sex)
          if(input_list$data$FishLenComps_LikeType[f] == 3) {
            map_FishLen_corr_pars[r,s,f,1] <- counter_fishlen_corr
            counter_fishlen_corr <- counter_fishlen_corr + 1
          }
        }

        # joint by sex, split by region (ages)
        if(any(fishage_comp_type == 2) && input_list$data$FishAgeComps_LikeType[f] != 0 && s == 1) {
          map_FishAge_theta[r,1,f] <- counter_fishage
          counter_fishage <- counter_fishage + 1 # joint by sex, split by region

          # if logistic normal 1dar1 by age only (unique age correlations for age each region and sex)
          if(input_list$data$FishAgeComps_LikeType[f] == 3) {
            map_FishAge_corr_pars[r,1,f,1] <- counter_fishage_corr
            counter_fishage_corr <- counter_fishage_corr + 1
          }

          # if logistic normal, 2d, where 1dar1 by age, constant correlation by sex
          if(input_list$data$FishAgeComps_LikeType[f] == 4) {
            for(i in 1:2) {
              if(i == 2 && input_list$data$n_sexes == 1) next # skip if we only have 1 sex
              map_FishAge_corr_pars[r,1,f,i] <- counter_fishage_corr
              counter_fishage_corr <- counter_fishage_corr + 1
            } # end i
          } # end if
        }

        # joint by sex, split by region (lengths)
        if(any(fishlen_comp_type == 2) && input_list$data$FishLenComps_LikeType[f] != 0 && s == 1) {
          map_FishLen_theta[r,1,f] <- counter_fishlen
          counter_fishlen <- counter_fishlen + 1 # joint by sex, split by region

          # if logistic normal 1dar1 by length only (unique length correlations for length each region and sex)
          if(input_list$data$FishLenComps_LikeType[f] == 3) {
            map_FishLen_corr_pars[r,1,f,1] <- counter_fishlen_corr
            counter_fishlen_corr <- counter_fishlen_corr + 1
          }

          # if logistic normal, 2d, where 1dar1 by length, constant correlation by sex
          if(input_list$data$FishLenComps_LikeType[f] == 4) {
            for(i in 1:2) {
              if(i == 2 && input_list$data$n_sexes == 1) next # skip if we only have 1 sex
              map_FishLen_corr_pars[r,1,f,i] <- counter_fishlen_corr
              counter_fishlen_corr <- counter_fishlen_corr + 1
            } # end i
          } # end if
        }

      } # end s loop

    } # end r loop
  } # end f loop

  # Input into mapping list
  input_list$map$ln_FishAge_theta <- factor(map_FishAge_theta)
  input_list$map$ln_FishLen_theta <- factor(map_FishLen_theta)
  input_list$map$ln_FishAge_theta_agg <- factor(map_FishAge_theta_agg)
  input_list$map$ln_FishLen_theta_agg <- factor(map_FishLen_theta_agg)
  input_list$map$FishAge_corr_pars_agg <- factor(map_FishAge_corr_pars_agg)
  input_list$map$FishAge_corr_pars <- factor(map_FishAge_corr_pars)
  input_list$map$FishLen_corr_pars_agg <- factor(map_FishLen_corr_pars_agg)
  input_list$map$FishLen_corr_pars <- factor(map_FishLen_corr_pars)

  # Print all messages if verbose is TRUE
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}

#' Setup fishery selectivity and catchability specifications
#'
#' @param input_list List containing a data list, parameter list, and map list
#' @param cont_tv_fish_sel Character vector specifying the form of continuous time-varying selectivity for each fishery fleet.
#' The vector must be length \code{n_fish_fleets}, and each element must follow the structure:
#' \code{"<time variation type>_Fleet_<fleet number>"}.
#'
#' Valid time variation types include:
#' \itemize{
#'   \item \code{"none"}: No continuous time variation.
#'   \item \code{"iid"}: Independent and identically distributed deviations across years.
#'   \item \code{"rw"}: Random walk in time.
#'   \item \code{"3dmarg"}: 3D marginal time-varying selectivity.
#'   \item \code{"3dcond"}: 3D conditional time-varying selectivity.
#'   \item \code{"2dar1"}: Two-dimensional AR1 process.
#' }
#'
#' For example:
#' \itemize{
#'   \item \code{"iid_Fleet_1"} applies an iid time-varying structure to Fleet 1.
#'   \item \code{"none_Fleet_2"} means no time variation is used for Fleet 2.
#' }
#'
#' \strong{Note:} If time-block-based selectivity (via \code{fish_sel_blocks}) is specified for a fleet, then its corresponding entry here must be \code{"none_Fleet_<fleet number>"}.
#' @param fish_sel_blocks Character vector specifying the fishery selectivity blocks for each region and fleet.
#' Each element must follow the structure: `"Block_<block number>_Year_<start>-<end>_Fleet_<fleet number>"` or `"none_Fleet_<fleet number>"`.
#' This allows users to define time-varying selectivity blocks for specific fleets within a region.
#'
#' For example:
#' \itemize{
#'   \item \code{"Block_1_Year_1-35_Fleet_1"} defines selectivity block 1 for Fleet 1 covering years 1 through 35.
#'   \item \code{"Block_2_Year_36-56_Fleet_1"} defines block 2 for Fleet 1 for years 36 to 56.
#'   \item \code{"Block_3_Year_57-terminal_Fleet_1"} assigns block 3 from year 57 through the terminal year for Fleet 1.
#'   \item \code{"none_Fleet_2"} indicates that no fishery selectivity blocks are used for Fleet 2.
#' }
#'
#' The blocks must be non-overlapping and sequential in time within each fleet, and block numbers must be unique within each fleet.
#' @param fish_sel_model Character vector specifying the fishery selectivity model for each fleet.
#' The vector must be length \code{n_fish_fleets}, and each element must follow the structure:
#' \code{"<selectivity model>_Fleet_<fleet number>"}.
#'
#' Available selectivity model types include:
#' \itemize{
#'   \item \code{"logist1"}: Logistic function with parameters \code{a50} and \code{k}.
#'   \item \code{"logist2"}: Logistic function with parameters \code{a50} and \code{a95}.
#'   \item \code{"gamma"}: Dome-shaped gamma function with parameters \code{amax} and \code{delta}.
#'   \item \code{"exponential"}: Exponential function with a power parameter.
#' }
#'
#' For example:
#' \itemize{
#'   \item \code{"logist1_Fleet_1"} uses the logistic (a50, k) model for Fleet 1.
#'   \item \code{"gamma_Fleet_2"} uses the gamma dome-shaped model for Fleet 2.
#' }
#'
#' The models are applied by region and year as defined in the overall model array structure
#' (\code{n_regions} x \code{n_years} x \code{n_fish_fleets}), though this vector defines only the functional form for each fleet.
#'
#' For mathematical definitions and implementation details of each selectivity form, refer to the model equations vignette.
#' @param fish_q_blocks Character vector specifying fishery catchability (q) blocks for each fleet.
#' Each element must follow the structure: \code{"Block_<block number>_Year_<start>-<end>_Fleet_<fleet number>"}
#' or \code{"none_Fleet_<fleet number>"}.
#'
#' This allows users to define time-varying catchability blocks independently of selectivity blocks.
#' The blocks must be non-overlapping and sequential in time within each fleet.
#'
#' For example:
#' \itemize{
#'   \item \code{"Block_1_Year_1-35_Fleet_1"} assigns block 1 to Fleet 1 for years 1–35.
#'   \item \code{"Block_2_Year_36-56_Fleet_1"} continues with block 2 for years 36–56.
#'   \item \code{"Block_3_Year_57-terminal_Fleet_1"} assigns block 3 from year 57 to the terminal year for Fleet 1.
#'   \item \code{"none_Fleet_2"} indicates no catchability blocks are used for Fleet 2.
#' }
#'
#' Internally, these specifications are converted to a \code{[n_regions, n_years, n_fish_fleets]} array,
#' where each block is mapped to the appropriate years and fleets.
#' @param fishsel_pe_pars_spec Character string specifying how process error parameters for fishery selectivity
#' are estimated across regions and sexes. This is only relevant if \code{cont_tv_fish_sel} is not set to \code{"none"};
#' otherwise, all process error parameters are treated as fixed.
#'
#' Available options include:
#' \itemize{
#'   \item \code{"est_all"}: Estimates separate process error parameters for each region and sex.
#'   \item \code{"est_shared_r"}: Shares process error parameters across regions (sex-specific parameters are still estimated).
#'   \item \code{"est_shared_s"}: Shares process error parameters across sexes (region-specific parameters are still estimated).
#'   \item \code{"est_shared_r_s"}: Shares process error parameters across both regions and sexes, estimating a single set of parameters.
#'   \item \code{"fix"} or \code{"none"}: Does not estimate process error parameters; all are treated as fixed.
#' }
#' @param fish_fixed_sel_pars_spec Character string specifying the structure for estimating
#' fixed-effect parameters of the fishery selectivity model (e.g., a50, k, amax).
#' This controls whether selectivity parameters are estimated separately or shared across regions and sexes.
#'
#' Available options include:
#' \itemize{
#'   \item \code{"est_all"}: Estimates separate fixed-effect selectivity parameters for each region and sex.
#'   \item \code{"est_shared_r"}: Shares parameters across regions (sex-specific parameters are still estimated).
#'   \item \code{"est_shared_s"}: Shares parameters across sexes (region-specific parameters are still estimated).
#'   \item \code{"est_shared_r_s"}: Shares parameters across both regions and sexes, estimating a single set of fixed-effect parameters.
#' }
#' @param fish_q_spec Character string specifying the structure of fishery catchability (\code{q}) estimation
#' across regions. This controls whether separate or shared parameters are used.
#'
#' Available options include:
#' \itemize{
#'   \item \code{"est_all"}: Estimates separate catchability parameters for each region.
#'   \item \code{"est_shared_r"}: Estimates a single catchability parameter shared across all regions.
#' }
#' @param fish_sel_devs_spec Character string specifying the structure of process error deviations
#' in time-varying fishery selectivity dimensioned by the number of fishery fleets. This determines how deviations are estimated across regions and sexes.
#'
#' Available options include:
#' \itemize{
#'   \item \code{"est_all"}: Estimates a separate deviation time series for each region and sex.
#'   \item \code{"est_shared_r"}: Shares deviations across regions (sex-specific deviations are still estimated).
#'   \item \code{"est_shared_s"}: Shares deviations across sexes (region-specific deviations are still estimated).
#'   \item \code{"est_shared_r_s"}: Shares deviations across both regions and sexes, estimating a single deviation time series.
#' }
#'
#' This argument is only used when a continuous time-varying selectivity form is specified (e.g., via \code{cont_tv_fish_sel}).
#' @param corr_opt_semipar Character string specifying which correlation structures to suppress
#'   when using semi-parametric time-varying selectivity models. Only used if \code{cont_tv_sel}
#'   is set to one of \code{"3dmarg"}, \code{"3dcond"}, or \code{"2dar1"}.
#'
#'   This option allows users to turn off estimation of specific correlation components in the
#'   time-varying selectivity model. This can improve stability or enforce assumptions about
#'   independence in the temporal or age structure.
#'
#'   Available options:
#'   \itemize{
#'     \item \code{"corr_zero_y"}: Sets year (temporal) correlations to 0.
#'     \item \code{"corr_zero_a"}: Sets age correlations to 0.
#'     \item \code{"corr_zero_y_a"}: Sets both year and age correlations to 0.
#'     \item \code{"corr_zero_c"}: Sets cohort correlations to 0. Only valid for \code{cont_tv_sel} = \code{"3dmarg"} or \code{"3dcond"}.
#'     \item \code{"corr_zero_y_c"}: Sets year and cohort correlations to 0. Only valid for \code{cont_tv_sel} = \code{"3dmarg"} or \code{"3dcond"}.
#'     \item \code{"corr_zero_a_c"}: Sets age and cohort correlations to 0. Only valid for \code{cont_tv_sel} = \code{"3dmarg"} or \code{"3dcond"}.
#'     \item \code{"corr_zero_y_a_c"}: Sets all correlations (year, age, and cohort) to 0.
#'       Only valid for \code{cont_tv_sel} = \code{"3dmarg"} or \code{"3dcond"}; equivalent to an iid structure.
#'   }
#'
#' These correlation-suppression flags are ignored when \code{cont_tv_sel} is set to any other value.
#' @param Use_fish_q_prior Integer specifying whether to use fishery q prior or not (0 dont use) (1 use)
#' @param fish_q_prior Fishery q priors in normal space, dimensioned by region, block,  fishery fleet, and 2 (mean, and sd in the 4 dimension of array)
#' @param ... Additional arguments specifying starting values for fishery selectivity and catchability parameters (fishsel_pe_pars, ln_fishsel_devs, ln_fish_fixed_sel_pars, ln_fish_q)
#'
#' @export Setup_Mod_Fishsel_and_Q
#'
#' @importFrom stringr str_detect
Setup_Mod_Fishsel_and_Q <- function(input_list,
                                    cont_tv_fish_sel,
                                    fish_sel_blocks,
                                    fish_sel_model,
                                    Use_fish_q_prior = 0,
                                    fish_q_prior = NA,
                                    fish_q_blocks,
                                    fishsel_pe_pars_spec = NULL,
                                    fish_fixed_sel_pars_spec = NULL,
                                    fish_q_spec = NULL,
                                    fish_sel_devs_spec = NULL,
                                    corr_opt_semipar = NULL,
                                    ...
                                    ) {

  messages_list <<- character(0) # string to attach to for printing messages

  if(!is.null(fishsel_pe_pars_spec)) if(length(fishsel_pe_pars_spec) != input_list$data$n_fish_fleets) stop("fishsel_pe_pars_spec is not length n_fish_fleets")
  if(!is.null(fish_sel_devs_spec)) if(length(fish_sel_devs_spec) != input_list$data$n_fish_fleets) stop("fish_sel_devs_spec is not length n_fish_fleets")
  if(!is.null(corr_opt_semipar)) if(length(corr_opt_semipar) != input_list$data$n_fish_fleets) stop("corr_opt_semipar is not length n_fish_fleets")
  if(!Use_fish_q_prior %in% c(0,1)) stop("Values for Use_fish_q_prior are not valid. They are == 0 (don't use prior), or == 1 (use prior)")
  collect_message("Fishery Catchability priors are: ", ifelse(Use_fish_q_prior == 0, "Not Used", "Used"))

  # define for continuous time-varying selectivity
  cont_tv_fish_sel_mat <- array(NA, dim = c(input_list$data$n_regions, input_list$data$n_fish_fleets))
  cont_tv_map <- data.frame(type = c("none", "iid", "rw", "3dmarg", "3dcond", "2dar1"),
                            num = c(0,1,2,3,4,5)) # set up values we map to

  for(i in 1:length(cont_tv_fish_sel)) {
    # Extract out components from list
    tmp <- cont_tv_fish_sel[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    cont_tv_type <- tmp_vec[1] # get continuous selex type
    fleet <- as.numeric(tmp_vec[3]) # extract fleet index
    if(!fleet %in% c(1:input_list$data$n_fish_fleets)) stop("Invalid fleet specified for cont_tv_fish_sel This needs to be specified as timevarytype_Fleet_x")
    if(!cont_tv_type %in% c(cont_tv_map$type)) stop("cont_tv_fish_sel is not correctly specified. This needs to be one of these: none, iid, rw, 3dmarg, 3dcond, 2dar1 (the timevarytypes) and specified as timevarytype_Fleet_x")
    cont_tv_fish_sel_mat[,fleet] <- cont_tv_map$num[which(cont_tv_map$type == cont_tv_type)]
    collect_message("Continuous fishery time-varying selectivity specified as: ", cont_tv_type, " for fishery fleet ", fleet)
  }

  # define fishery time blocks
  fish_sel_blocks_arr = array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets))
  for(i in 1:length(fish_sel_blocks)) {
    # Extract out components from list
    tmp <- fish_sel_blocks[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    if(!tmp_vec[1] %in% c("none", "Block")) stop("Fishery Selectivity Blocks not correctly specified. This should be either none_Fleet_x or Block_x_Year_x-y_Fleet_x")
    # extract out fleets if constant
    if(tmp_vec[1] == "none") {
      fleet <- as.numeric(tmp_vec[3]) # get fleet number
      fish_sel_blocks_arr[,,fleet] <- 1 # input only 1 fishery time block
    }
    if(tmp_vec[1] == "Block") {
      block_val <- as.numeric(tmp_vec[2]) # get block value

      # get year ranges
      if(!str_detect(tmp, "terminal")) { # if not terminal year
        year_range <- as.numeric(unlist(strsplit(tmp_vec[4], "-")))
        years <- year_range[1]:year_range[2] # get sequence of years
      } else { # if terminal year
        year_range <- unlist(strsplit(tmp_vec[4], '-'))[1] # get year range
        years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
      }

      fleet <- as.numeric(tmp_vec[6]) # extract fleet index
      fish_sel_blocks_arr[,years,fleet] <- block_val
    }
  }

  if(any(is.na(fish_sel_blocks_arr))) stop("Fishery Selectivtiy Blocks are returning an NA. Did you forget to specify the year range of fish_sel_blocks?")

  for(f in 1:input_list$data$n_fish_fleets) collect_message(paste("Fishery Selectivity Time Blocks for fishery", f, "is specified at:", length(unique(fish_sel_blocks_arr[,,f]))))

  # Setup fishery selectivity models (functional forms)
  sel_map <- data.frame(sel = c('logist1', "gamma", "exponential", "logist2"), num = c(0,1,2,3)) # set up values we can map to
  fish_sel_model_arr = array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets))
  for(i in 1:length(fish_sel_model)) {
    # Extract out components from list
    tmp <- fish_sel_model[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    sel_type <- tmp_vec[1] # get selectivity type
    fleet <- as.numeric(tmp_vec[3]) # get fleet number
    if(!sel_type %in% c(sel_map$sel)) stop("fish_sel_model is not correctly specified. This needs to be one of these: logist1, gamma, exponential, logist2 (the seltypes) and specified as seltype_Fleet_x")
    if(!fleet %in% c(1:input_list$data$n_fish_fleets)) stop("Invalid fleet specified for fish_sel_model This needs to be specified as seltype_Fleet_x")
    fish_sel_model_arr[,,fleet] <- sel_map$num[which(sel_map$sel == sel_type)]
    collect_message("Fishery selectivity functional form specified as:", sel_type, " for fishery fleet ", f)
  } # end i loop

  # setup fishery catchability blocks
  fish_q_blocks_arr = array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets))
  for(i in 1:length(fish_q_blocks)) {
    # Extract out components from list
    tmp <- fish_q_blocks[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    if(!tmp_vec[1] %in% c("none", "Block")) stop("Fishery Catchability Blocks not correctly specified. This should be either none_Fleet_x or Block_x_Year_x-y_Fleet_x")
    # extract out fleets if constant
    if(tmp_vec[1] == "none") {
      fleet <- as.numeric(tmp_vec[3]) # get fleet number
      fish_q_blocks_arr[,,fleet] <- 1 # input only 1 fishery catchability time block
    }
    if(tmp_vec[1] == "Block") {
      block_val <- as.numeric(tmp_vec[2]) # get block value

      # get year ranges
      if(!str_detect(tmp, "terminal")) { # if not terminal year
        year_range <- as.numeric(unlist(strsplit(tmp_vec[4], "-")))
        years <- year_range[1]:year_range[2] # get sequence of years
      } else { # if terminal year
        year_range <- unlist(strsplit(tmp_vec[4], '-'))[1] # get year range
        years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
      }

      fleet <- as.numeric(tmp_vec[6]) # get fleet number
      fish_q_blocks_arr[,years,fleet] <- block_val # input catchability time block
    }
  }

  if(any(is.na(fish_q_blocks))) stop("Fishery Catchability Blocks are returning an NA. Did you forget to specify the year range of fish_q_blocks?")

  for(f in 1:input_list$data$n_fish_fleets) collect_message(paste("Fishery Catchability Time Blocks for fishery", f, "is specified at:", length(unique(fish_q_blocks_arr[,,f]))))

  # Setup data input list
  input_list$data$cont_tv_fish_sel <- cont_tv_fish_sel_mat
  input_list$data$fish_sel_blocks <- fish_sel_blocks_arr
  input_list$data$fish_sel_model <- fish_sel_model_arr
  input_list$data$fish_q_blocks <- fish_q_blocks_arr
  input_list$data$fish_q_prior <- fish_q_prior
  input_list$data$Use_fish_q_prior <- Use_fish_q_prior

  # Set up parameter inputs
  starting_values <- list(...)

  # Fishery selectivity fixed effects
  # Figure out number of selectivity parameters for a given functional form
  unique_fishsel_vals <- unique(as.vector(input_list$data$fish_sel_model))
  sel_pars_vec <- vector() # create empty vector to populate

  for(i in 1:length(unique_fishsel_vals)) {
    if(unique_fishsel_vals[i] %in% c(2)) sel_pars_vec[i] <- 1
    if(unique_fishsel_vals[i] %in% c(0,1,3)) sel_pars_vec[i] <- 2
  } # end i loop

  max_fishsel_blks <- max(apply(input_list$data$fish_sel_blocks, c(1,3), FUN = function(x) length(unique(x)))) # figure out maximum number of fishery selectivity blocks for a given reigon and fleet
  max_fishsel_pars <- max(sel_pars_vec) # maximum number of selectivity parameters across all forms
  if("ln_fish_fixed_sel_pars" %in% names(starting_values)) input_list$par$ln_fish_fixed_sel_pars <- starting_values$ln_fish_fixed_sel_pars
  else input_list$par$ln_fish_fixed_sel_pars <- array(0, dim = c(input_list$data$n_regions, max_fishsel_pars, max_fishsel_blks, input_list$data$n_sexes, input_list$data$n_fish_fleets))

  # Fishery catchability
  max_fishq_blks <- max(apply(input_list$data$fish_q_blocks, c(1,3), FUN = function(x) length(unique(x)))) # figure out maximum number of fishery catchability blocks for a given reigon and fleet
  if("ln_fish_q" %in% names(starting_values)) input_list$par$ln_fish_q <- starting_values$ln_fish_q
  else input_list$par$ln_fish_q <- array(0, dim = c(input_list$data$n_regions, max_fishq_blks, input_list$data$n_fish_fleets))

  # Fishery selectivity process error parameters
  if("fishsel_pe_pars" %in% names(starting_values)) input_list$par$fishsel_pe_pars <- starting_values$fishsel_pe_pars
  else input_list$par$fishsel_pe_pars <- array(0, dim = c(input_list$data$n_regions, max(max_fishsel_pars, 4), input_list$data$n_sexes, input_list$data$n_fish_fleets)) # dimensioned 4 as the max number of pars for process errors (e.g., sigmas), and then just map off if not using

  # Fishery selectivity deviations
  if("ln_fishsel_devs" %in% names(starting_values)) input_list$par$ln_fishsel_devs <- starting_values$ln_fishsel_devs
  else input_list$par$ln_fishsel_devs <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), length(input_list$data$ages), input_list$data$n_sexes, input_list$data$n_fish_fleets))

  # Setup mapping list
  # Initialize counter and mapping array for fixed effects fishery selectivity
  fish_fixed_sel_pars_counter <- 1
  map_fish_fixed_sel_pars <- input_list$par$ln_fish_fixed_sel_pars
  map_fish_fixed_sel_pars[] <- NA

  for(f in 1:input_list$data$n_fish_fleets) {
    for(r in 1:input_list$data$n_regions) {

      if(!fish_fixed_sel_pars_spec[f] %in% c("est_all", "est_shared_r", "est_shared_r_s", "fix", "est_shared_s"))
        stop("fish_fixed_sel_pars_spec not correctly specfied. Should be one of these: est_all, est_shared_r, est_shared_r_s, est_shared_s, fix")

      # Figure out max number of selectivity parameters for a given region and fleet
      if(unique(input_list$data$fish_sel_model[r,,f]) %in% 2) max_sel_pars <- 1 # exponential
      if(unique(input_list$data$fish_sel_model[r,,f]) %in% c(0,1,3)) max_sel_pars <- 2 # logistic a50, k, gamma, and logistic a50, a95

      # Extract number of fishery selectivity blocks
      fishsel_blocks_tmp <- unique(as.vector(input_list$data$fish_sel_blocks[r,,f]))

      for(s in 1:input_list$data$n_sexes) {
        for(b in 1:length(fishsel_blocks_tmp)) {
          for(i in 1:max_sel_pars) {

            # Estimate all selectivity fixed effects parameters within the constraints of the defined blocks
            if(fish_fixed_sel_pars_spec[f] == "est_all") {
              map_fish_fixed_sel_pars[r,i,b,s,f] <- fish_fixed_sel_pars_counter
              fish_fixed_sel_pars_counter <- fish_fixed_sel_pars_counter + 1
            } # end if

            # Estimating process error parameters shared across regions (but unique for each sex, fleet, parameter)
            if(fish_fixed_sel_pars_spec[f] == 'est_shared_r' && r == 1) {
              for(rr in 1:input_list$data$n_regions) {
                if(fishsel_blocks_tmp[b] %in% input_list$data$fish_sel_blocks[rr,,f]) {
                  map_fish_fixed_sel_pars[rr, i, b, s, f] <- fish_fixed_sel_pars_counter
                } # end if
              } # end rr loop
              fish_fixed_sel_pars_counter <- fish_fixed_sel_pars_counter + 1
            } # end if

            # Estimating process error parameters shared across sexes (but unique for each region, fleet, parameter)
            if(fish_fixed_sel_pars_spec[f] == 'est_shared_s' && s == 1) {
              for(ss in 1:input_list$data$n_sexes) {
                map_fish_fixed_sel_pars[r, i, b, ss, f] <- fish_fixed_sel_pars_counter
              } # end ss loop
              fish_fixed_sel_pars_counter <- fish_fixed_sel_pars_counter + 1
            } # end if

            # Estimating process error parameters shared across regions and sexes (but unique for each fleet, parameter)
            if(fish_fixed_sel_pars_spec[f] == 'est_shared_r_s' && r == 1 && s == 1) {
              for(rr in 1:input_list$data$n_regions) {
                for(ss in 1:input_list$data$n_sexes) {
                  if(fishsel_blocks_tmp[b] %in% input_list$data$fish_sel_blocks[rr,,f]) {
                    map_fish_fixed_sel_pars[rr, i, b, ss, f] <- fish_fixed_sel_pars_counter
                  } # end if
                } # end ss loop
              } #end rr loop
              fish_fixed_sel_pars_counter <- fish_fixed_sel_pars_counter + 1
            } # end if

          } # end i loop
        } # end b loop
      } # end s loop
    } # end r loop

    # fix all parameters
    if(fish_fixed_sel_pars_spec[f] == "fix") map_fish_fixed_sel_pars[,,,,f] <- NA
    collect_message("fish_fixed_sel_pars_spec is specified as: ", fish_fixed_sel_pars_spec[f], " for fishery fleet ", f)
  } # end f loop

  # Initialize counter and mapping array for fishery catchability
  fish_q_counter <- 1
  map_fish_q <- input_list$par$ln_fish_q
  map_fish_q[] <- NA

  for(f in 1:input_list$data$n_fish_fleets) {

    if(!is.null(fish_q_spec)) {
      if(!fish_q_spec[f] %in% c("est_all", "est_shared_r", "fix"))
        stop("fish_q_spec not correctly specfied. Should be one of these: est_all, est_shared_r, fix")
    }

    for(r in 1:input_list$data$n_regions) {
      if(sum(input_list$data$UseFishIdx[r,,f]) == 0) {
        map_fish_q[r,,f] <- NA # fix parameters if we are not using fishery indices for these fleets and regions
      } else {
        # Extract number of fishery catchability blocks
        fishq_blocks_tmp <- unique(as.vector(input_list$data$fish_q_blocks[r,,f]))
        for(b in 1:length(fishq_blocks_tmp)) {

          # Estimate for all regions
          if(fish_q_spec[f] == 'est_all') {
           map_fish_q[r,b,f] <- fish_q_counter
           fish_q_counter <- fish_q_counter + 1
          } # end if

          # Estimate but share q across regions
          if(fish_q_spec[f] == 'est_shared_r' && r == 1) {
            for(rr in 1:input_list$data$n_regions) {
              if(fishq_blocks_tmp[b] %in% input_list$data$fish_q_blocks[rr,,f]) {
                map_fish_q[rr, b, f] <- fish_q_counter
              } # end if
            } # end rr loop
            fish_q_counter <- fish_q_counter + 1
          } # end if

        } # end b loop
        # fix all parameters
        if(fish_q_spec[f] == 'fix') map_fish_q[,,f] <- NA
      } # end else loop
    } # end r loop
    collect_message("fish_q_spec is specified as: ", fish_q_spec[f], " for fishery fleet ", f)
  } # end f loop

  # Initialize counter and mapping array for fishery process errors
  fishsel_pe_pars_counter <- 1 # initalize counter
  map_fishsel_pe_pars <- input_list$par$fishsel_pe_pars # initalize array
  map_fishsel_pe_pars[] <- NA

  # Fishery process error parameters
  for(f in 1:input_list$data$n_fish_fleets) {

    if(!is.null(fishsel_pe_pars_spec)) {
      if(!fishsel_pe_pars_spec[f] %in% c("fix", "none", "est_all", "est_shared_r", "est_shared_s", "est_shared_r_s"))
        stop("fishsel_pe_pars_spec not correctly specfied. Should be one of these: fix, none, est_all, est_shared_r, est_shared_r_s, est_shared_s")
    }

    for(r in 1:input_list$data$n_regions) {

      # Figure out max number of selectivity parameters for a given region and fleet
      if(unique(input_list$data$fish_sel_model[r,,f]) %in% 2) max_sel_pars <- 1 # exponential
      if(unique(input_list$data$fish_sel_model[r,,f]) %in% c(0,1,3)) max_sel_pars <- 2 # logistic a50, k, gamma, and logistic a50, a95

      # if no time-variation, then fix all parameters for this fleet
      if(input_list$data$cont_tv_fish_sel[r,f] == 0) {
        map_fishsel_pe_pars[r,,,f] <- NA
      } else { # if we have time-variation
        for(s in 1:input_list$data$n_sexes) {
          # If iid time-variation or random walk for this fleet
          if(input_list$data$cont_tv_fish_sel[r,f] %in% c(1,2)) {
            for(i in 1:max_sel_pars) {
              # either fixing parameters or not used for a given fleet
              if(fishsel_pe_pars_spec[f] %in% c("none", "fix")) map_fishsel_pe_pars[r,i,s,f] <- NA
              # Estimating all parameters separately (unique for each region, sex, fleet, parameter)
              if(fishsel_pe_pars_spec[f] == "est_all") {
                map_fishsel_pe_pars[r,i,s,f] <- fishsel_pe_pars_counter
                fishsel_pe_pars_counter <- fishsel_pe_pars_counter + 1
              } # end est_all
              # Estimating process error parameters shared across regions (but unique for each sex, fleet, parameter)
              if(fishsel_pe_pars_spec[f] == 'est_shared_r' && r == 1) {
                map_fishsel_pe_pars[,i,s,f] <- fishsel_pe_pars_counter
                fishsel_pe_pars_counter <- fishsel_pe_pars_counter + 1
              }
              # Estimating process error parameters shared across sexes (but unique for each region, fleet, parameter)
              if(fishsel_pe_pars_spec[f] == 'est_shared_s' && s == 1) {
                map_fishsel_pe_pars[r,i,,f] <- fishsel_pe_pars_counter
                fishsel_pe_pars_counter <- fishsel_pe_pars_counter + 1
              }
              # Estimating process error parameters shared across regions and sexes (but unique for each fleet, parameter)
              if(fishsel_pe_pars_spec[f] == 'est_shared_r_s' && r == 1 && s == 1) {
                map_fishsel_pe_pars[,i,,f] <- fishsel_pe_pars_counter
                fishsel_pe_pars_counter <- fishsel_pe_pars_counter + 1
              }
            } # end i loop
          } # end iid or random walk variation

          # If 3d gmrf or 2dar1
          if(input_list$data$cont_tv_fish_sel[r,f] %in% c(3,4,5)) {

            # Set up indexing to loop through
            if(input_list$data$cont_tv_fish_sel[r,f] %in% c(3,4)) idx = 1:4 # 3dgmrf (1 = pcorr_age, 2 = pcorr_year, 3= pcorr_cohort, 4 = log_sigma)
            if(input_list$data$cont_tv_fish_sel[r,f] %in% c(5)) idx = c(1,2,4) # 2dar1 (1 = pcorr_age, 2 = pcorr_year, 4 = log_sigma)

            for(i in idx) {

              # either fixing parameters or not used for a given fleet
              if(fishsel_pe_pars_spec[f] %in% c("none", "fix")) map_fishsel_pe_pars[r,i,s,f] <- NA

              # Estimating all process error parameters
              if(fishsel_pe_pars_spec[f] == "est_all") {
                map_fishsel_pe_pars[r,i,s,f] <- fishsel_pe_pars_counter
                fishsel_pe_pars_counter <- fishsel_pe_pars_counter + 1
              } # end est_all
              # Estimating process error parameters shared across regions (but unique for each sex, fleet, parameter)
              if(fishsel_pe_pars_spec[f] == 'est_shared_r' && r == 1) {
                map_fishsel_pe_pars[,i,s,f] <- fishsel_pe_pars_counter
                fishsel_pe_pars_counter <- fishsel_pe_pars_counter + 1
              }
              # Estimating process error parameters shared across sexes (but unique for each region, fleet, parameter)
              if(fishsel_pe_pars_spec[f] == 'est_shared_s' && s == 1) {
                map_fishsel_pe_pars[r,i,,f] <- fishsel_pe_pars_counter
                fishsel_pe_pars_counter <- fishsel_pe_pars_counter + 1
              }
              # Estimating process error parameters shared across regions and sexes (but unique for each fleet, parameter)
              if(fishsel_pe_pars_spec[f] == 'est_shared_r_s' && r == 1 && s == 1) {
                map_fishsel_pe_pars[,i,,f] <- fishsel_pe_pars_counter
                fishsel_pe_pars_counter <- fishsel_pe_pars_counter + 1
              }
            } # end i loop

            # Options to set correaltions to 0 for 3dgmrf
            if(!is.null(corr_opt_semipar)) {

              if(!corr_opt_semipar[f] %in% c(NA, "corr_zero_y", "corr_zero_a", "corr_zero_y_a", "corr_zero_c", "corr_zero_y_c", "corr_zero_a_c", "corr_zero_y_a_c"))
                stop("corr_opt_semipar not correctly specfied. Should be one of these: corr_zero_y, corr_zero_a, corr_zero_y_a, corr_zero_c, corr_zero_y_c, corr_zero_a_c, corr_zero_y_a_c, NA")

              # if either 2dar1 or 3dgmrf
              if(input_list$data$cont_tv_fish_sel[r,f] %in% c(3,4,5)) {
                if(corr_opt_semipar[f] == "corr_zero_y") map_fishsel_pe_pars[,2,,1] <- NA
                if(corr_opt_semipar[f] == "corr_zero_a") map_fishsel_pe_pars[,1,,1] <- NA
                if(corr_opt_semipar[f] == "corr_zero_y_a") map_fishsel_pe_pars[,1:2,,1] <- NA
              }

              # if 3dgmrf only
              if(input_list$data$cont_tv_fish_sel[r,f] %in% c(3,4)) {
                if(corr_opt_semipar[f] == "corr_zero_c") map_fishsel_pe_pars[,3,,1] <- NA
                if(corr_opt_semipar[f] == "corr_zero_y_c") map_fishsel_pe_pars[,2:3,,1] <- NA
                if(corr_opt_semipar[f] == "corr_zero_a_c") map_fishsel_pe_pars[,c(1,3),,1] <- NA
                if(corr_opt_semipar[f] == "corr_zero_y_a_c") map_fishsel_pe_pars[,1:3,,1] <- NA
              }

              # Reset numbering for mapping off correlation parameters for clarity
              non_na_positions <- which(!is.na(map_fishsel_pe_pars))
              map_fishsel_pe_pars[non_na_positions] <- seq_along(non_na_positions)
              collect_message("corr_opt_semipar is specified as: ", corr_opt_semipar[f], "for fishery fleet", f)

            }

          } # end if 3d gmrf marginal or conditional variance

          # fix all parameters
          if(fishsel_pe_pars_spec[f] == "fix") map_fishsel_pe_pars[r,,s,f] <- NA

        } # end s loop
      } # end else
    } # end r loop

    if(!is.null(fishsel_pe_pars_spec)) collect_message("fishsel_pe_pars_spec is specified as: ", fishsel_pe_pars_spec[f], "for fishery fleet", f)

  } # end f loop

  # Initialize counter and mapping array for fishery selectivity deviations
  fishsel_devs_counter <- 1
  map_fishsel_devs <- input_list$par$ln_fishsel_devs
  map_fishsel_devs[] <- NA

  for(r in 1:input_list$data$n_regions) {
    for(f in 1:input_list$data$n_fish_fleets) {

      if(!is.null(fish_sel_devs_spec)) {
        if(!fish_sel_devs_spec[f] %in% c("none", "est_all", "est_shared_r", "est_shared_s", "est_shared_r_s"))
          stop("fish_sel_devs_spec not correctly specfied. Should be one of these: est_all, est_shared_r, est_shared_r_s, none")
      }

    # Figure out max number of selectivity parameters for a given region and fleet
    if(unique(input_list$data$fish_sel_model[r,,f]) %in% 2) max_sel_pars <- 1 # exponential
    if(unique(input_list$data$fish_sel_model[r,,f]) %in% c(0,1,3)) max_sel_pars <- 2 # logistic a50, k, gamma, and logistic a50, a95

      for(s in 1:input_list$data$n_sexes) {
        for(y in 1:length(input_list$data$years)) {

          # if no time-variation, then fix all parameters for this fleet
          if(input_list$data$cont_tv_fish_sel[r,f] == 0) {
            map_fishsel_devs[r,y,,s,f] <- NA
          } else {
            # If iid or random walk time-variation for this fleet
            if(input_list$data$cont_tv_fish_sel[r,f] %in% c(1,2) && y >= min(which(input_list$data$UseCatch[,,f] == 1))) {
              for(i in 1:max_sel_pars) {
                # Estimating all selectivity deviations across regions, sexes, fleets, and parameter
                if(fish_sel_devs_spec[f] == 'est_all') {
                  map_fishsel_devs[r,y,i,s,f] <- fishsel_devs_counter
                  fishsel_devs_counter <- fishsel_devs_counter + 1
                }
                # Estimating selectivity deviations across sexes, fleets, and parameters, but shared across regions
                if(fish_sel_devs_spec[f] == 'est_shared_r' && r == 1) {
                  map_fishsel_devs[,y,i,s,f] <- fishsel_devs_counter
                  fishsel_devs_counter <- fishsel_devs_counter + 1
                }
                # Estimating selectivity deviations across regions, fleets, and parameters, but shared across sexes
                if(fish_sel_devs_spec[f] == 'est_shared_s' && s == 1) {
                  map_fishsel_devs[r,y,i,,f] <- fishsel_devs_counter
                  fishsel_devs_counter <- fishsel_devs_counter + 1
                }
                # Estimating selectivity deviations across fleets, and parameters, but shared across sexes and regions
                if(fish_sel_devs_spec[f] == 'est_shared_r_s' && r == 1 && s == 1) {
                  map_fishsel_devs[,y,i,,f] <- fishsel_devs_counter
                  fishsel_devs_counter <- fishsel_devs_counter + 1
                }
              } # end i loop
            } # end iid or random walk variation

            # If 3d gmrf for this fleet
            if(input_list$data$cont_tv_fish_sel[r,f] %in% c(3,4,5) && y >= min(which(input_list$data$UseCatch[,,f] == 1))) {
              for(i in 1:length(input_list$data$ages)) {
                # Estimating all selectivity deviations across regions, years and ages (also cohorts baked in year x age)
                if(fish_sel_devs_spec[f] == 'est_all') {
                  map_fishsel_devs[r,y,i,s,f] <- fishsel_devs_counter
                  fishsel_devs_counter <- fishsel_devs_counter + 1
                }
                # Estimating all selectivity deviations across years and ages (also cohorts baked in year x age), but shared across regions
                if(fish_sel_devs_spec[f] == 'est_shared_r' && r == 1) {
                  map_fishsel_devs[,y,i,s,f] <- fishsel_devs_counter
                  fishsel_devs_counter <- fishsel_devs_counter + 1
                }
                # Estimating all selectivity deviations across years and ages (also cohorts baked in year x age), but shared across sexes
                if(fish_sel_devs_spec[f] == 'est_shared_s' && s == 1) {
                  map_fishsel_devs[r,y,i,,f] <- fishsel_devs_counter
                  fishsel_devs_counter <- fishsel_devs_counter + 1
                }

                # Estimating all selectivity deviations across years and ages (also cohorts baked in year x age), but shared across sexes and regions
                if(fish_sel_devs_spec[f] == 'est_shared_r_s' && s == 1 && r == 1) {
                  map_fishsel_devs[,y,i,,f] <- fishsel_devs_counter
                  fishsel_devs_counter <- fishsel_devs_counter + 1
                }

              } # end i loop
            } # end 3d gmrf

          } # end else
        } # end y loop
      } # end s loop

    if(!is.null(fish_sel_devs_spec)) collect_message("fish_sel_devs_spec is specified as: ", fish_sel_devs_spec[f], "for fishery fleet", f, "and region ", r)

    } # end f loop
  } # end r loop

  map_fishsel_devs <<- map_fishsel_devs

  # input into mapping list
  input_list$map$fishsel_pe_pars <- factor(map_fishsel_pe_pars)
  input_list$map$ln_fish_fixed_sel_pars <- factor(map_fish_fixed_sel_pars)
  input_list$map$ln_fish_q <- factor(map_fish_q)
  input_list$map$ln_fishsel_devs <- factor(map_fishsel_devs)
  input_list$data$map_ln_fishsel_devs <- array(as.numeric(input_list$map$ln_fishsel_devs), dim = dim(input_list$par$ln_fishsel_devs))
  input_list$data$map_fish_q <- map_fish_q

  # Checking whether fishery q dimensions are correct
  if(Use_fish_q_prior == 1) if(sum(dim(fish_q_prior) == c(dim(map_fish_q), 2)) != 4) stop("Fishery catchability dimensions are not correct. Should be n_regions, max n_blocks, n_fish_fleets, and 2 (where 2 represents the 2 prior parameters - the mean and sd). You can input an NA if not availiable for certain regions or fleets.")

  # Print all messages if verbose is TRUE
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}


