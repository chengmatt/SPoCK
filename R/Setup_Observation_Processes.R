#' Set up simulation observation processes
#'
#' @param Comp_Structure Composition Structure (spltR_spltS, spltR_jntS)
#' @param Comp_Srv_Like Survey Composition Likelihoods (Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal)
#' @param Srv_Like_Pars Parameters for Survey Composition Likelihoods
#' @param Comp_Fish_Like Fishery Composition Likelihoods (Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal)
#' @param Tag_Like Tag Likelihoods (Poisson, NegBin, Multinomial_Release, Multinomial_Recapture)
#' @param Tag_Like_Pars Tag Likelihood Parameters
#' @param base_ISS_FishAge Base Input sample size for fishery ages
#' @param sim_list Simulation list
#' @param base_ISS_SrvAge Input sample size for survey ages
#' @param ISS_FishAge_Pattern Input sample size pattern for fishery ages. Either "constant" or follows an f pattern ("F_pattern")
#' @param SrvAgeTheta Survey Age Composition Overdispersion Parameter
#' @param FishAgeTheta Fishery Age Composition Overdispersion Parameter
#'
#' @export Setup_Sim_Observation_Proc
#'
Setup_Sim_Observation_Proc <- function(Comp_Structure,
                                       Comp_Srv_Like,
                                       base_ISS_SrvAge,
                                       Srv_Like_Pars,
                                       Comp_Fish_Like,
                                       base_ISS_FishAge,
                                       ISS_FishAge_Pattern,
                                       Tag_Like,
                                       Tag_Like_Pars,
                                       SrvAgeTheta,
                                       FishAgeTheta,
                                       sim_list
                                       ) {
  # Create empty containers
  comp_srv_like <- vector()
  comp_fish_like <- vector()

  # Set up composition structure
  if(Comp_Structure == "spltR_spltS") sim_list$comp_strc <- 0
  if(Comp_Structure == "spltR_jntS") sim_list$comp_strc <- 1

  # Set up survey composition likelihoods
  for(sf in 1:sim_list$n_srv_fleets) {
    if(Comp_Srv_Like[sf] == "Multinomial") comp_srv_like <- c(comp_srv_like, 0)
    if(Comp_Srv_Like[sf] == "Dirichlet-Multinomial") comp_srv_like <- c(comp_srv_like, 1)
  } # end sf loop


  # Set up fishery composition likelihoods
  for(f in 1:sim_list$n_fish_fleets) {
    if(Comp_Fish_Like[f] == "Multinomial") comp_fish_like <- c(comp_fish_like, 0)
    if(Comp_Fish_Like[f] == "Dirichlet-Multinomial") comp_fish_like <- c(comp_fish_like, 1)
    if(Comp_Fish_Like[f] == "iid-Logistic-Normal") comp_fish_like <- c(comp_fish_like, 2)
    if(Comp_Fish_Like[f] == "3d-Logistic-Normal") comp_fish_like <- c(comp_fish_like, 5)
  } # end f loop

  # Set up tagging likelihoods
  if(Tag_Like == "Poisson") sim_list$tag_like <- 0
  if(Tag_Like == "NegBin") sim_list$tag_like <- 1
  if(Tag_Like == "Multinomial_Release") sim_list$tag_like <- 2
  if(Tag_Like == "Multinomial_Recapture") sim_list$tag_like <- 3

  # setup variance and overdispersion parameters for composition and tagging data
  if(ISS_FishAge_Pattern == 'constant') ISS_FishAge <- array(base_ISS_FishAge, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_fish_fleets, sim_list$n_sims))

  # output into list
  sim_list$Tag_Like_Pars <- Tag_Like_Pars
  sim_list$Srv_Like_Pars <- Srv_Like_Pars
  sim_list$comp_srv_like <- comp_srv_like
  sim_list$comp_fish_like <- comp_fish_like
  sim_list$ISS_FishAge <- ISS_FishAge
  sim_list$ISS_SrvAge <- base_ISS_SrvAge
  sim_list$FishAgeTheta <- FishAgeTheta
  sim_list$SrvAgeTheta <- SrvAgeTheta

  return(sim_list)
}

#' Set up SPoRC model weighting
#'
#' @param input_list List containing data, parameter, and map lists.
#' @param likelihoods Numeric flag indicating likelihood implementation to use:
#'   \itemize{
#'     \item 0 for ADMB likelihoods
#'     \item 1 for TMB likelihoods
#'   }
#' @param Wt_Catch Either a numeric scalar (lambda) applied to the overall catch dataset or an array of lambdas (i.e., weights can change by year and fleet) dimensioned by n_regions, n_years, n_fish_fleets.
#' @param Wt_FishIdx Either a numeric scalar (lambda) applied to the overall fishery index dataset  or an array of lambdas (i.e., weights can change by year and fleet) dimensioned by n_regions, n_years, n_fish_fleets.
#' @param Wt_SrvIdx Either a numeric scalar (lambda) applied to the overall survey index dataset or an array of lambdas (i.e., weights can change by year and fleet) dimensioned by n_regions, n_years, n_srv_fleets.
#' @param Wt_Rec Numeric weight (lambda) applied to the recruitment penalty.
#' @param Wt_F Numeric weight (lambda) applied to the fishing mortality penalty.
#' @param Wt_FishAgeComps Numeric weight (lambda) applied to fishery age composition data.
#' @param Wt_SrvAgeComps Numeric weight (lambda) applied to survey age composition data.
#' @param Wt_FishLenComps Numeric weight (lambda) applied to fishery length composition data.
#' @param Wt_SrvLenComps Numeric weight (lambda) applied to survey length composition data.
#' @param sablefish_ADMB Numeric flag to mimic calculations for the sablefish ADMB model:
#'   \itemize{
#'     \item 1 to mimic sablefish ADMB calculations
#'     \item 0 otherwise
#'   }
#' @param Wt_Tagging Numeric weight (lambda) applied to tagging data.
#'
#' @export Setup_Mod_Weighting
Setup_Mod_Weighting <- function(input_list,
                                sablefish_ADMB,
                                likelihoods,
                                Wt_Catch = 1,
                                Wt_FishIdx = 1,
                                Wt_SrvIdx = 1,
                                Wt_Rec = 1,
                                Wt_F = 1,
                                Wt_Tagging = 1,
                                Wt_FishAgeComps,
                                Wt_SrvAgeComps,
                                Wt_FishLenComps,
                                Wt_SrvLenComps
                                ) {

  messages_list <<- character(0) # string to attach to for printing messages

  if(!likelihoods %in% c(0,1)) stop("likelihoods are not correctly specified. Should be either 0 (ADMB) or 1 (TMB)")
  else collect_message("Using ", ifelse(likelihoods == 0, 'ADMB', 'TMB'), " likelihoods")

  input_list$data$sablefish_ADMB <- sablefish_ADMB
  input_list$data$likelihoods <- likelihoods
  input_list$data$Wt_Catch <- Wt_Catch
  input_list$data$Wt_FishIdx <- Wt_FishIdx
  input_list$data$Wt_SrvIdx <- Wt_SrvIdx
  input_list$data$Wt_Rec <- Wt_Rec
  input_list$data$Wt_F <- Wt_F
  input_list$data$Wt_FishAgeComps<- Wt_FishAgeComps
  input_list$data$Wt_SrvAgeComps<- Wt_SrvAgeComps
  input_list$data$Wt_FishLenComps<- Wt_FishLenComps
  input_list$data$Wt_SrvLenComps<- Wt_SrvLenComps
  input_list$data$Wt_Tagging <- Wt_Tagging

  # Print all messages if verbose is TRUE
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}
