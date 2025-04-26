#' Set up simulation observaiton processes
#'
#' @param n_fish_fleets Number of fishery fleets
#' @param n_srv_fleets Number of survey fleets
#' @param Comp_Structure Composition Structure (SpltSex_SpltRegion, JntSex_SpltRegion, JntSex_JntRegion)
#' @param Comp_Srv_Like Survey Composition Likelihoods (Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal)
#' @param Srv_Like_Pars Parameters for Survey Composition Likelihoods
#' @param Comp_Fish_Like Fishery Composition Likelihoods (Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal)
#' @param Fish_Like_Pars Parameters for Fishery Composition Likelihoods
#' @param Tag_Like Tag Likelihoods (Poisson, NegBin, Multinomial_Release, Multinomial_Recapture)
#' @param Tag_Like_Pars Tag Likelihood Parameters
#'
#' @export Setup_Sim_Observation_Proc
#'
Setup_Sim_Observation_Proc <- function(n_fish_fleets,
                                       n_srv_fleets,
                                       Comp_Structure,
                                       Comp_Srv_Like,
                                       Srv_Like_Pars,
                                       Comp_Fish_Like,
                                       Fish_Like_Pars,
                                       Tag_Like,
                                       Tag_Like_Pars) {
  # Create empty containers
  comp_srv_like <- vector()
  comp_fish_like <- vector()

  # Set up composition structure
  if(Comp_Structure == "SpltSex_SpltRegion") comp_strc <<- 0
  if(Comp_Structure == "JntSex_SpltRegion") comp_strc <<- 1
  if(Comp_Structure == "JntSex_JntRegion") comp_strc <<- 2

  # Set up survey composition likelihoods
  for(sf in 1:n_srv_fleets) {
    if(Comp_Srv_Like[sf] == "Multinomial") comp_srv_like <- c(comp_srv_like, 0)
    if(Comp_Srv_Like[sf] == "Dirichlet-Multinomial") comp_srv_like <- c(comp_srv_like, 1)
  } # end sf loop


  # Set up fishery composition likelihoods
  for(f in 1:n_fish_fleets) {
    if(Comp_Fish_Like[f] == "Multinomial") comp_fish_like <- c(comp_fish_like, 0)
    if(Comp_Fish_Like[f] == "Dirichlet-Multinomial") comp_fish_like <- c(comp_fish_like, 1)
    if(Comp_Fish_Like[f] == "iid-Logistic-Normal") comp_fish_like <- c(comp_fish_like, 2)
    if(Comp_Fish_Like[f] == "3d-Logistic-Normal") comp_fish_like <- c(comp_fish_like, 5)
  } # end f loop

  # Set up tagging likelihoods
  if(Tag_Like == "Poisson") tag_like <<- 0
  if(Tag_Like == "NegBin") tag_like <<- 1
  if(Tag_Like == "Multinomial_Release") tag_like <<- 2
  if(Tag_Like == "Multinomial_Recapture") tag_like <<- 3

  # output into globals
  Tag_Like_Pars <<- Tag_Like_Pars
  Fish_Like_Pars <<- Fish_Like_Pars
  Srv_Like_Pars <<- Srv_Like_Pars
  comp_srv_like <<- comp_srv_like
  comp_fish_like <<- comp_fish_like
}

#' Title
#'
#' @param input_list List containing a data list, parameter list, and map list
#' @param likelihoods Numeric indicating whether to use ADMB likelihoods (0) or TMB likelihoods (1)
#' @param Wt_Catch Weight (lambda) applied to the overall catch dataset
#' @param Wt_FishIdx Weight (lambda) applied to the overall fishery index dataset
#' @param Wt_SrvIdx Weight (lambda) applied to the overall survey index dataset
#' @param Wt_Rec Weight (lambda) applied to recruitment penalty
#' @param Wt_F Weight (lambda) applied to fishing mortality penalty
#' @param Wt_FishAgeComps Weight (lambda) applied to fishery age compositions
#' @param Wt_SrvAgeComps Weight (lambda) applied to survey age compositions
#' @param Wt_FishLenComps Weight (lambda) applied to fishery length compositions
#' @param Wt_SrvLenComps Weight (lambda) applied to survey length compositions
#' @param sablefish_ADMB Numeric indicating whether to mimic calculations for the sablefish ADMB model (1) or not (0)
#'
#' @export Setup_Mod_Weighting
#'
Setup_Mod_Weighting <- function(input_list,
                                sablefish_ADMB,
                                likelihoods,
                                Wt_Catch = 1,
                                Wt_FishIdx = 1,
                                Wt_SrvIdx = 1,
                                Wt_Rec = 1,
                                Wt_F = 1,
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

  # Print all messages if verbose is TRUE
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}
