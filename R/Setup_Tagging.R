#' Set up simulated tagging dynamics
#'
#' @param n_tags Number of tags to release in a given year
#' @param max_liberty Maximum liberty to track cohorts for
#' @param tag_years Years to release tags for
#' @param t_tagging Time of tagging (e.g., start year == 0, mid year == 0.5)
#' @param base_Tag_Reporting Base tag reporting rate by region
#' @param Tag_Reporting_pattern Tag reporting pattern. Options include: constant
#' @param Tag_Ind_Mort Initial tag induced mortality
#' @param Tag_Shed Chronic tag shedding rate
#' @param sim_list Simulation list
#'
#' @export Setup_Sim_Tagging
#'
Setup_Sim_Tagging <- function(
                              n_tags,
                              max_liberty,
                              tag_years,
                              t_tagging,
                              base_Tag_Reporting,
                              Tag_Reporting_pattern,
                              Tag_Ind_Mort,
                              Tag_Shed,
                              sim_list
                              ) {

  # Output variables into list
  sim_list$n_tags <- n_tags
  sim_list$max_liberty <- max_liberty
  sim_list$tag_years <- tag_years
  sim_list$n_tag_yrs <- length(tag_years)
  sim_list$t_tagging <- t_tagging
  sim_list$Tag_Ind_Mort <- Tag_Ind_Mort
  sim_list$Tag_Shed <- Tag_Shed
  sim_list$n_tag_rel_events <- sim_list$n_tag_yrs * sim_list$n_regions # number of tag release events - tag years x tag region (tag cohorts)
  sim_list$tag_rel_indicator <- expand.grid(regions = 1:sim_list$n_regions, tag_yrs = tag_years) # get tag release indicator (by tag years and regions = a tag cohort)

  # Containers
  Tag_Reporting <- array(0, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims)) # populated later on
  sim_list$Tag_Fish <- array(0, dim = c(sim_list$n_tag_rel_events, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims)) # number of tagged fish
  sim_list$Tag_Avail <- array(0, dim = c(sim_list$max_liberty + 1, sim_list$n_tag_rel_events, sim_list$n_regions, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims)) # tags availiable for recapture every year
  sim_list$Pred_Tag_Recap <- array(0, dim = c(sim_list$max_liberty, sim_list$n_tag_rel_events, sim_list$n_regions, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims)) # predicted tag recaptures
  sim_list$Obs_Tag_Recap <- array(0, dim = c(sim_list$max_liberty, sim_list$n_tag_rel_events, sim_list$n_regions, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims)) # observed tag recaptures

  for(sim in 1:sim_list$n_sims) {
    for(r in 1:sim_list$n_regions) {
      for(y in 1:sim_list$n_yrs) {
        if(Tag_Reporting_pattern == "constant") Tag_Reporting[r,y,sim] <- base_Tag_Reporting[r]
      } # end y loop
    } # end r loop
  } # end sim loop

  sim_list$Tag_Reporting <- Tag_Reporting # output this into list

  return(sim_list)
}

#' Setup tagging processes and parameters
#'
#' @param input_list List containing a data list, parameter list, and map list
#' @param UseTagging Numeric (0 or 1) indicating whether to use tagging data (1) or not (0)
#' @param tag_release_indicator Matrix [n_tag_cohorts x 2], where columns are release region and release year
#' @param max_tag_liberty Maximum number of years to track a tagged cohort
#' @param Tagged_Fish Array [n_tag_cohorts x n_ages x n_sexes] describing tagged fish releases
#' @param Obs_Tag_Recap Array [max_tag_liberty x n_tag_cohorts x n_regions x n_ages x n_sexes] observed tag recaptures
#' @param Tag_LikeType Character string specifying tag likelihood type. One of:
#'   \itemize{
#'     \item \code{"Poisson"}
#'     \item \code{"NegBin"}
#'     \item \code{"Multinomial_Release"}
#'     \item \code{"Multinomial_Recapture"}
#'     \item \code{"Dirichlet-Multinomial_Release"}
#'     \item \code{"Dirichlet-Multinomial_Recapture"}
#'   }
#'   Example: \code{Tag_LikeType = "NegBin"}
#' @param mixing_period Numeric indicating minimum years post-release to include in fitting
#' @param t_tagging Fractional year when tagging occurs (e.g., 0.5 for mid-year)
#' @param tag_selex Character string specifying tag recovery selectivity. One of:
#'   \itemize{
#'     \item \code{"Uniform_DomFleet"}
#'     \item \code{"SexAgg_DomFleet"}
#'     \item \code{"SexSp_DomFleet"}
#'     \item \code{"Uniform_AllFleet"}
#'     \item \code{"SexAgg_AllFleet"}
#'     \item \code{"SexSp_AllFleet"}
#'   }
#'   Example: \code{tag_selex = "SexSp_AllFleet"}
#' @param tag_natmort Character string specifying tag natural mortality parameterization. One of:
#'   \itemize{
#'     \item \code{"AgeAgg_SexAgg"}
#'     \item \code{"AgeSp_SexAgg"}
#'     \item \code{"AgeAgg_SexSp"}
#'     \item \code{"AgeSp_SexSp"}
#'   }
#'   Example: \code{tag_natmort = "AgeSp_SexSp"}
#' @param Use_TagRep_Prior Numeric (0 or 1) whether to use tag reporting rate prior
#' @param TagRep_PriorType Numeric indicating prior type for tag reporting:
#'   \itemize{
#'     \item 0: Symmetric beta
#'     \item 1: Regular beta
#'   }
#' @param TagRep_mu Numeric mean for tag reporting prior (normal space); \code{NA} if symmetric beta is used
#' @param TagRep_sd Numeric standard deviation for tag reporting prior (normal space)
#' @param move_age_tag_pool List or character specifying pooling of tagging data by age groups. Examples:
#'   \itemize{
#'     \item \code{list(c(1:6), c(7:15), c(16:30))} pools ages 1–6, 7–15, 16–30
#'     \item \code{"all"} pools all ages together as one group (internally converted to \code{list(1:n_ages)})
#'     \item \code{list(1, 2, 3, ..., n_ages)} fits each age separately
#'   }
#' @param move_sex_tag_pool List or character specifying pooling of tagging data by sex groups. Examples:
#'   \itemize{
#'     \item \code{list(1:2)} pools sexes together
#'     \item \code{"all"} pools all sexes together (internally converted to \code{list(1:n_sexes)})
#'     \item \code{list(1, 2)} fits each sex separately
#'   }
#' @param Init_Tag_Mort_spec Character string \code{"fix"} or \code{"est"} specifying if initial tag mortality is fixed or estimated
#' @param Tag_Shed_spec Character string \code{"fix"} or \code{"est"} specifying if chronic tag shedding is fixed or estimated
#' @param Tag_Reporting_blocks Character vector specifying blocks of years and regions for tag reporting rates. Format examples:
#'   \itemize{
#'     \item \code{"Block_1_Year_1-15_Region_1"}
#'     \item \code{"Block_2_Year_16-terminal_Region_2"}
#'     \item \code{"none_Region_3"} (means no block, constant for that region)
#'   }
#' @param TagRep_spec Character string specifying tag reporting rate estimation scheme:
#'   \itemize{
#'     \item \code{"est_all"} estimates rates for all blocks and regions independently
#'     \item \code{"est_shared_r"} estimates rates shared across regions but varying by block
#'     \item \code{"fix"} fixes all reporting rates (no estimation)
#'   }
#' @param ... Additional starting values for tagging parameters such as \code{ln_Init_Tag_Mort}, \code{ln_Tag_Shed}, \code{ln_tag_theta}, \code{Tag_Reporting_Pars}
#'
#' @export Setup_Mod_Tagging
Setup_Mod_Tagging <- function(input_list,
                              UseTagging = 0,
                              tag_release_indicator = NULL,
                              max_tag_liberty = 0,
                              Tagged_Fish = NA,
                              Obs_Tag_Recap = NA,
                              Tag_LikeType = NA,
                              mixing_period = 1,
                              t_tagging = 0,
                              tag_selex = NA,
                              tag_natmort = NA,
                              Use_TagRep_Prior = 0,
                              TagRep_PriorType = NA,
                              TagRep_mu = NA,
                              TagRep_sd = NA,
                              move_age_tag_pool = NA,
                              move_sex_tag_pool = NA,
                              Init_Tag_Mort_spec = NULL,
                              Tag_Shed_spec = NULL,
                              TagRep_spec = NULL,
                              Tag_Reporting_blocks = NULL,
                              ...
                              ) {

  messages_list <<- character(0) # string to attach to for printing messages

  # Setup tagging likelihood
  tag_like_map <- data.frame(type = c("Poisson", "NegBin", "Multinomial_Release", "Multinomial_Recapture", "Dirichlet-Multinomial_Release", "Dirichlet-Multinomial_Recapture"), num = c(0,1,2,3,4,5))
  if(is.na(Tag_LikeType)) Tag_LikeType_vals <- 999
  else {
    if(!Tag_LikeType %in% c(tag_like_map$type)) stop("Tag Likelihood not correctly specified. Should be one of these: Poisson, NegBin, Multinomial_Release, Multinomial_Recapture, Dirichlet-Multinomial_Release, Dirichlet-Multinomial_Recapture")
    Tag_LikeType_vals <- tag_like_map$num[tag_like_map$type == Tag_LikeType]
    collect_message("Tag Likelihood specified as: ", Tag_LikeType)
  }

  # Setup tagging selectivity
  tag_selex_map <- data.frame(type = c("Uniform_DomFleet", "SexAgg_DomFleet", "SexSp_DomFleet",
                                       "Uniform_AllFleet", "SexAgg_AllFleet", "SexSp_AllFleet"), num = c(0,1,2,3,4,5))
  if(is.na(Tag_LikeType)) tag_selex_vals <- 999
  else {
    if(!tag_selex %in% c(tag_selex_map$type)) stop("Tag Selectivity not correctly specified. Should be one of these: Uniform_DomFleet, SexAgg_DomFleet, SexSp_DomFleet, Uniform_AllFleet, SexAgg_AllFleet, SexSp_AllFleet")
    tag_selex_vals <- tag_selex_map$num[tag_selex_map$type == tag_selex]
    collect_message("Tag Selectivity specified as: ", tag_selex)
  }

  # setup tagging natural mortality
  tag_natmort_map <- data.frame(type = c("AgeAgg_SexAgg", "AgeSp_SexAgg", "AgeAgg_SexSp", "AgeSp_SexSp"), num = c(0,1,2,3))
  if(is.na(Tag_LikeType)) tag_natmort_vals <- 999
  else {
    if(!tag_natmort %in% c(tag_natmort_map$type)) stop("Tag Natural Mortality not correctly specified. Should be one of these: AgeAgg_SexAgg, AgeSp_SexAgg, AgeAgg_SexSp, AgeSp_SexSp")
    tag_natmort_vals <- tag_natmort_map$num[tag_natmort_map$type == tag_natmort]
     collect_message("Tag Natural Mortality specified as: ", tag_natmort)
  }

  # If movement is pooled either across sexes or ages
  if(is.character(move_age_tag_pool)){
    if(move_age_tag_pool == "all") move_age_tag_pool_vals = list(input_list$data$ages)
  } else move_age_tag_pool_vals = move_age_tag_pool

  if(is.character(move_sex_tag_pool)){
    if(move_sex_tag_pool == "all") move_sex_tag_pool_vals = list(1:input_list$data$n_sexes)
  } else move_sex_tag_pool_vals = move_sex_tag_pool

  collect_message("Tagging data are fit to by pooling across ", length(move_age_tag_pool_vals), " age groups")
  collect_message("Tagging data are fit to by pooling across ", length(move_sex_tag_pool_vals), " sex groups")

  # setup tag reporting rates
  Tag_Reporting_blocks_mat = array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years)))
  if(!is.null(Tag_Reporting_blocks)) {
    for(i in 1:length(Tag_Reporting_blocks)) {
      # Extract out components from list
      tmp <- Tag_Reporting_blocks[i]
      tmp_vec <- unlist(strsplit(tmp, "_"))
      if(!tmp_vec[1] %in% c("none", "Block")) stop("Tag Reporting Blocks not correctly specified. This should be either none_Region_x or Block_x_Year_x-y_Region_x")
      # extract out fleets if constant
      if(tmp_vec[1] == "none") {
        region <- as.numeric(tmp_vec[3]) # get region index
        Tag_Reporting_blocks_mat[region,] <- 1 # input tag reporting time block
      }
      if(tmp_vec[1] == "Block") {
        block_val <- as.numeric(tmp_vec[2]) # get block value
        region <- as.numeric(tmp_vec[6]) # get region value

        # get year ranges
        if(!str_detect(tmp, "terminal")) { # if not terminal year
          year_range <- as.numeric(unlist(strsplit(tmp_vec[4], "-")))
          years <- year_range[1]:year_range[2] # get sequence of years
        } else { # if terminal year
          year_range <- unlist(strsplit(tmp_vec[4], '-'))[1] # get year range
          years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
        }

        Tag_Reporting_blocks_mat[region,years] <- block_val # input tag reporting time block
      }
    }
  }

   for(r in 1:input_list$data$n_regions) collect_message("Tag Reporting estimated with ", length(unique(Tag_Reporting_blocks_mat[r,])), " block for region ", r)

  # Input data list
  input_list$data$UseTagging <- UseTagging
  input_list$data$tag_release_indicator <- tag_release_indicator
  if(UseTagging == 0) input_list$data$n_tag_cohorts <- 0
  if(UseTagging == 1) input_list$data$n_tag_cohorts <- nrow(tag_release_indicator)
  input_list$data$max_tag_liberty <- max_tag_liberty
  input_list$data$Tagged_Fish <- Tagged_Fish
  input_list$data$Obs_Tag_Recap <- Obs_Tag_Recap
  input_list$data$Tag_LikeType <- Tag_LikeType_vals
  input_list$data$mixing_period <- mixing_period
  input_list$data$t_tagging <- t_tagging
  input_list$data$tag_selex <- tag_selex_vals
  input_list$data$tag_natmort <- tag_natmort_vals
  input_list$data$Use_TagRep_Prior <- Use_TagRep_Prior
  input_list$data$TagRep_PriorType <- TagRep_PriorType
  input_list$data$TagRep_mu <- TagRep_mu
  input_list$data$TagRep_sd <- TagRep_sd
  input_list$data$move_age_tag_pool <- move_age_tag_pool_vals
  input_list$data$move_sex_tag_pool <- move_sex_tag_pool_vals
  input_list$data$Tag_Reporting_blocks <- Tag_Reporting_blocks_mat

  # Input parameter list
  starting_values <- list(...)

  # Initial tag induced mortality
  if("ln_Init_Tag_Mort" %in% names(starting_values)) input_list$par$ln_Init_Tag_Mort <- starting_values$ln_Init_Tag_Mort
  else input_list$par$ln_Init_Tag_Mort <- log(1e-50)

  # Chronic tag shedding
  if("ln_Tag_Shed" %in% names(starting_values)) input_list$par$ln_Tag_Shed <- starting_values$ln_Tag_Shed
  else input_list$par$ln_Tag_Shed <- log(1e-50)

  # tag overdispersion parameter
  if("ln_tag_theta" %in% names(starting_values)) input_list$par$ln_tag_theta <- starting_values$ln_tag_theta
  else input_list$par$ln_tag_theta <- 0

  # tag reporting rate parameters
  max_tagrep_blks <- max(apply(input_list$data$Tag_Reporting_blocks, 1, FUN = function(x) length(unique(x)))) # figure out maximum number of tag reporting rate blocks for each region
  if("Tag_Reporting_Pars" %in% names(starting_values)) input_list$par$Tag_Reporting_Pars <- starting_values$Tag_Reporting_Pars
  else input_list$par$Tag_Reporting_Pars <- array(0, dim = c(input_list$data$n_regions, max_tagrep_blks)) # specified at 0.5 in inverse logit space

  # Setup mapping list
  # If not using tagging data
  if(UseTagging == 0) {
    input_list$map$ln_Init_Tag_Mort <- factor(NA) # initial tag mortality
    input_list$map$ln_Tag_Shed <- factor(NA) # chronic tag shedding
    input_list$map$ln_tag_theta <- factor(NA) # tag overdispersion
    input_list$map$Tag_Reporting_Pars <- factor(rep(NA, length(input_list$par$Tag_Reporting_Pars))) # tag reporting rates
    input_list$data$map_Tag_Reporting_Pars = array(as.numeric(input_list$map$Tag_Reporting_Pars), dim = dim(input_list$par$Tag_Reporting_Pars))
  }

  # if using tagging data
  if(UseTagging == 1) {

    if(!Init_Tag_Mort_spec %in% c("fix", "est")) stop("Init_Tag_Mort_spec is incorrectly specified. Should be one of these: fix, est")
    if(!Tag_Shed_spec %in% c("fix", "est")) stop("Tag_Shed_spec is incorrectly specified. Should be one of these: fix, est")

    collect_message("Initial Tag Mortality is specified as: ", Init_Tag_Mort_spec)
    collect_message("Chronic Tag Shedding is specified as: ", Tag_Shed_spec)

    # Initial tag mortality
    if(Init_Tag_Mort_spec == "fix") input_list$map$ln_Init_Tag_Mort <- factor(NA)
    if(Init_Tag_Mort_spec == "est") input_list$map$ln_Init_Tag_Mort <- factor(1)

    # Tag Shedding
    if(Tag_Shed_spec == "fix" || UseTagging == 0) input_list$map$ln_Tag_Shed <- factor(NA)
    if(Tag_Shed_spec == "est") input_list$map$ln_Tag_Shed <- factor(1)

    # Tag Overdispersion
    if(input_list$data$Tag_LikeType %in% c(0,2,3)) input_list$map$ln_tag_theta <- factor(NA)
    if(input_list$data$Tag_LikeType %in% c(1,4,5)) input_list$map$ln_tag_theta <- factor(1)

    # Tag Reporting Rates
    # Initialize arrays and counters
    map_TagRep <- input_list$par$Tag_Reporting_Pars
    map_TagRep[] <- NA
    tagrep_counter <- 1

    # If this is a poisson, negative binomial, multinomial or dirichlet-multinomial release conditioned
    if(input_list$data$Tag_LikeType %in% c(0,1,2,4)) {
      for(r in 1:input_list$data$n_regions) {
        if(!is.null(TagRep_spec)) if(!TagRep_spec %in% c("est_all", "est_shared_r", "fix")) stop("Tag Reporting Specificaiton is not correctly specified. Needs to be fix, est_all, or est_shared_r")
        # Get number of tag reporting rate blocks
        tagrep_blocks_tmp <- unique(as.vector(input_list$data$Tag_Reporting_blocks[r,]))
        for(b in 1:length(tagrep_blocks_tmp)) {
          # Estimate for all regions
          if(TagRep_spec == 'est_all') {
            map_TagRep[r,b] <- tagrep_counter
            tagrep_counter <- tagrep_counter + 1
          }
          # Estimate but share tag reporting across regions
          if(TagRep_spec == 'est_shared_r' && r == 1) {
            for(rr in 1:input_list$data$n_regions) {
              # only assign if this value exists for this region
              if(tagrep_blocks_tmp[b] %in% input_list$data$Tag_Reporting_blocks[rr,]) {
                map_TagRep[rr, b] <- tagrep_counter
              } # end if
            } # end rr loop
            tagrep_counter <- tagrep_counter + 1
          }
        } # end b loop
      } # end r loop
      # if we want to fix
      if(TagRep_spec == 'fix') map_TagRep[] <- NA
      collect_message("Tag Reporting is specified as: ", TagRep_spec)
    } # end if

    # input tag reporting rates into mapping list
    input_list$map$Tag_Reporting_Pars <- factor(map_TagRep) # tag reporting rates
    input_list$data$map_Tag_Reporting_Pars <- array(as.numeric(input_list$map$Tag_Reporting_Pars), dim = dim(input_list$par$Tag_Reporting_Pars))

  } # end if for using tagging data

  # Print all messages if verbose is TRUE
  if(input_list$verbose && UseTagging == 1) for(msg in messages_list) message(msg)

  return(input_list)
}
