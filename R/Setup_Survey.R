#' Setup values for survey parameterization
#'
#' @param n_sims Number of simulations
#' @param n_yrs Number of years
#' @param n_regions Number of regions
#' @param n_ages Number of ages
#' @param n_sexes Number of sexes
#' @param n_srv_fleets Number of survey fleets
#' @param sigmaSrvIdx Survey index observation error, dimensioned by region and fleet
#' @param base_srv_q Base survey catchability value, dimensioned by region and fleet
#' @param srv_q_pattern Survey catchability pattern, dimensioned by region and fleet. Options include: constant
#' @param sel_model Survey selectivity model dimensioned by region and fleet. Options include: logistic
#' @param fixed_srv_sel_pars Fixed parameters of survey selectivity, dimensioned by region, sex, survey fleet, and the max number of parameters needed for
#' for a defined survey selectivity functional form out of all defined functional forms for the survey.
#'
#' @export Setup_Sim_Survey
#'
Setup_Sim_Survey <- function(n_sims = n_sims,
                             n_yrs = n_yrs,
                             n_regions = n_regions,
                             n_ages = n_ages,
                             n_sexes = n_sexes,
                             n_srv_fleets = n_srv_fleets,
                             sigmaSrvIdx,
                             base_srv_q,
                             srv_q_pattern,
                             sel_model,
                             fixed_srv_sel_pars
                             ) {

  # otuput survey sigma into globals
  sigmaSrvIdx <<- sigmaSrvIdx

  # create containers to loop through and populate
  srv_sel <- array(0, dim = c(n_yrs, n_regions, n_ages, n_sexes, n_srv_fleets, n_sims)) # survey selectivity
  srv_q <- array(1, dim = c(n_yrs, n_regions, n_srv_fleets, n_sims)) # survey catchability

  for(sim in 1:n_sims) {
    for(r in 1:n_regions) {
      for(y in 1:n_yrs) {
        for(sf in 1:n_srv_fleets) {

          # Survey catchability if constant
          if(srv_q_pattern[r,sf] == 'constant') srv_q[y,r,sf,sim] <- base_srv_q[r,sf]

          for(s in 1:n_sexes) {

            if(sel_model[r,sf] == 'logistic') {
              a50 <- fixed_srv_sel_pars[r,s,sf,1] # get a50
              k <- fixed_srv_sel_pars[r,s,sf,2] # get k
              srv_sel[y,r,,s,sf,sim] <- 1 / (1 + exp(-k * (1:n_ages - a50)))
            } # end if logistic

          }
        } # end f loop
      } # end y loop
    } # end r loop
  } # end sim loop

  # output into global environment
  srv_sel <<- srv_sel
  srv_q <<- srv_q


} # end function

#' Setup observed survey indices and composition data (age and length comps)
#'
#' @param input_list List containing a data list, parameter list, and map list
#' @param ObsSrvIdx Observed survey indices as an array dimensioned by n_regions, n_years, n_srv_fleets
#' @param ObsSrvIdx_SE Observed standard errors for survey indices as an array dimensioned by n_regions, n_years, n_srv_fleets
#' @param UseSrvIdx Indicator variables specifying whether or not to fit survey indices as an array imensioned by n_regions, n_years, n_srv_fleets, == 0 don't fit, == 1 fit
#' @param ObsSrvAgeComps Observed survey age compositions as an array dimensioned by n_regions, n_years, n_ages, n_sexes, n_srv_fleets (should be in numbers, although need not be whole numbers)
#' @param UseSrvAgeComps  Indicator variables specifying whether or not to fit survey ages as an array imensioned by n_regions, n_years, n_srv_fleets, == 0 don't fit, == 1 fit
#' @param ObsSrvLenComps Observed survey age compositions as an array dimensioned by n_regions, n_years, n_lens, n_sexes, n_srv_fleets (should be in numbers, although need not be whole numbers)
#' @param UseSrvLenComps  Indicator variables specifying whether or not to fit survey lengths as an array imensioned by n_regions, n_years, n_srv_fleets, == 0 don't fit, == 1 fit
#' @param SrvAgeComps_LikeType Vector of character strings dimensioned by n_srv_fleets, options include Multinomial, Dirichlet-Multinomial, and iid-Logistic-Normal. Can specify with "none" if no likelihoods are used.
#' @param SrvLenComps_LikeType Vector of character strings dimensioned by n_srv_fleets, options include Multinomial, Dirichlet-Multinomial, and iid-Logistic-Normal. Can specify with "none" if no likelihoods are used.
#' @param SrvAgeComps_Type Character vector defines the compositiont and that is structured as: the composition type (agg, spltRspltS, spltRjntS, jntRjntS), _, Year, _, year range you want the composition type to be in, _, Fleet, _, Fleet number. e.g., ("agg_Year_1-10_Fleet_1", "jntRjntS_Year_11-20_Fleet_1") species compositions to be aggregated in years 1 - 10 for fleet 1, joint by region and sex in years 11 - 20 for fleet 1.
#' @param SrvLenComps_Type Character vector defines the compositiont and that is structured as: the composition type (agg, spltRspltS, spltRjntS, jntRjntS, none), _, Year, _, year range you want the composition type to be in, _, Fleet, _, Fleet number. e.g., ("agg_Year_1-10_Fleet_1", "jntRjntS_Year_11-20_Fleet_1") species compositions to be aggregated in years 1 - 10 for fleet 1, joint by region and sex in years 11 - 20 for fleet 1.
#' @param SrvAge_comp_agg_type Vector dimensioned by n_srv_fleets specifying how to aggregate age composition data if SrvAgeComps_Type == 0, if comp_agg_type == 0 normalize, aggregate, ageing erorr then normalize, == 1 aggregate, normalize, then ageing error (for age compositions). Default is specify at NULL
#' @param SrvLen_comp_agg_type Vector dimensioned by n_srv_fleets specifying how to aggregate age composition data if SrvLenComps_Type == 1, if comp_agg_type == 0 length compositions are not normalized prior to application of size age transition, == 1 length compositions are normalized and then a size-age transition is applied (for length compositions). Default is specify at NULL
#' @param srv_idx_type Survey index type dimensioned by n_srv_fleets, character string for abundance (abd) or biomass (biom), "none" if not availabile or used
#' @param ... Additional arguments specifying starting values for survey age and length overdispersion parameters (ln_SrvAge_theta. ln_SrvLen_theta, ln_SrvAge_theta_agg, ln_SrvLen_theta_agg)
#' @param ISS_SrvAgeComps Input sample size for survey comps. Dimensioned by n_regions, n_years, n_sexes, n_fleets. Can be left as is if numbers are defined for the observed comps because the function sums them up interally to get an ISS. However, if observed comps are input as proportions, be sure to specify an ISS so the model doesn't incorrectly weight it as a value of ISS = 1.
#' @param ISS_SrvLenComps Input sample size for survey comps. Dimensioned by n_regions, n_years, n_sexes, n_fleets. Can be left as is if numbers are defined for the observed comps because the function sums them up interally to get an ISS. However, if observed comps are input as proportions, be sure to specify an ISS so the model doesn't incorrectly weight it as a value of ISS = 1.
#'
#' @export Setup_Mod_SrvIdx_and_Comps
#' @import stringr
Setup_Mod_SrvIdx_and_Comps <- function(input_list,
                                       ObsSrvIdx,
                                       ObsSrvIdx_SE,
                                       UseSrvIdx,
                                       srv_idx_type,
                                       ObsSrvAgeComps,
                                       UseSrvAgeComps,
                                       ObsSrvLenComps,
                                       UseSrvLenComps,
                                       ISS_SrvAgeComps = NULL,
                                       ISS_SrvLenComps  = NULL,
                                       SrvAgeComps_LikeType,
                                       SrvLenComps_LikeType,
                                       SrvAgeComps_Type,
                                       SrvLenComps_Type,
                                       SrvAge_comp_agg_type = NULL,
                                       SrvLen_comp_agg_type = NULL,
                                       ...
                                       ) {

  messages_list <<- character(0) # string to attach to for printing messages

  # Checking dimensions
  check_data_dimensions(ObsSrvIdx, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_srv_fleets = input_list$data$n_srv_fleets, what = 'ObsSrvIdx')
  check_data_dimensions(ObsSrvIdx_SE, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_srv_fleets = input_list$data$n_srv_fleets, what = 'ObsSrvIdx_SE')
  check_data_dimensions(UseSrvIdx, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_srv_fleets = input_list$data$n_srv_fleets, what = 'UseSrvIdx')
  check_data_dimensions(ObsSrvAgeComps, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_ages = length(input_list$data$ages), n_sexes = input_list$data$n_sexes, n_srv_fleets = input_list$data$n_srv_fleets, what = 'ObsSrvAgeComps')
  check_data_dimensions(UseSrvAgeComps, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_srv_fleets = input_list$data$n_srv_fleets, what = 'UseSrvAgeComps')
  check_data_dimensions(UseSrvLenComps, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_srv_fleets = input_list$data$n_srv_fleets, what = 'UseSrvLenComps')
  check_data_dimensions(ObsSrvLenComps, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_lens = length(input_list$data$lens), n_sexes = input_list$data$n_sexes, n_srv_fleets = input_list$data$n_srv_fleets, what = 'ObsSrvLenComps')
  if(!is.null(ISS_SrvAgeComps)) check_data_dimensions(ISS_SrvAgeComps, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_srv_fleets = input_list$data$n_srv_fleets, what = 'ISS_SrvAgeComps')
  if(!is.null(ISS_SrvLenComps)) check_data_dimensions(ISS_SrvLenComps, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_srv_fleets = input_list$data$n_srv_fleets, what = 'ISS_SrvLenComps')
  check_data_dimensions(srv_idx_type, n_srv_fleets = input_list$data$n_srv_fleets, what = 'srv_idx_type')
  check_data_dimensions(SrvAgeComps_LikeType, n_srv_fleets = input_list$data$n_srv_fleets, what = 'SrvAgeComps_LikeType')
  check_data_dimensions(SrvLenComps_LikeType, n_srv_fleets = input_list$data$n_srv_fleets, what = 'SrvLenComps_LikeType')

  # Checking specifications for index, age comps, length comps
  if(!all(srv_idx_type %in% c("biom", "abd", "none"))) stop("Invalid specification for srv_idx_type. Should be either abd, biom, or none")
  if(!all(SrvAgeComps_LikeType %in% c("none", "Multinomial", "Dirichlet-Multinomial", "iid-Logistic-Normal", "1d-Logistic-Normal", "2d-Logistic-Normal", "3d-Logistic-Normal")))
    stop("Invalid specification for SrvAgeComps_LikeType Should be either none, Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal, 1d-Logistic-Normal, 2d-Logistic-Normal, 3d-Logistic-Normal")
  if(!all(SrvLenComps_LikeType %in% c("none", "Multinomial", "Dirichlet-Multinomial", "iid-Logistic-Normal", "1d-Logistic-Normal", "2d-Logistic-Normal", "3d-Logistic-Normal")))
    stop("Invalid specification for SrvLenComps_LikeType Should be either none, Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal, 1d-Logistic-Normal, 2d-Logistic-Normal, 3d-Logistic-Normal")

  # Checking whether or not data are avalilable for all regions if jntR is specified
  if(!all(unique(as.vector(apply(UseSrvAgeComps, 2, sum))) %in% c(0, input_list$data$n_regions)) &&
     any(str_detect(SrvAgeComps_Type, "jntR"))) warning("Age composition data are not availiable in certain regions for a given year. It is not recommended to structure composition data as jntR because the bin sizes change, which could impact inference")
  if(!all(unique(as.vector(apply(UseSrvLenComps, 2, sum))) %in% c(0, input_list$data$n_regions)) &&
     any(str_detect(SrvLenComps_Type, "jntR"))) warning("Length composition data are not availiable in certain regions for a given year. It is not recommended to structure composition data as jntR because the bin sizes change, which could impact inference")

  # Specifying survey index type
  srv_idx_type_vals <- array(NA, dim = c(input_list$data$n_regions, input_list$data$n_srv_fleets))
  for(f in 1:ncol(srv_idx_type_vals)) {
    if(srv_idx_type[f] == 'biom') srv_idx_type_vals[,f] <- 1 # biomass
    if(srv_idx_type[f] == 'abd') srv_idx_type_vals[,f] <- 0 # abundance
    if(srv_idx_type[f] == 'none') srv_idx_type_vals[,f] <- 999 # none
    collect_message(paste("Survey Index", "for survey fleet", f, "specified as:" , srv_idx_type[f]))
  } # end f loop

  # Specifying survey age composition values
  comp_srvage_like_vals <- vector()
  for(f in 1:input_list$data$n_srv_fleets) {
    if(SrvAgeComps_LikeType[f] == 'none') comp_srvage_like_vals <- c(comp_srvage_like_vals, 999)
    if(SrvAgeComps_LikeType[f] == "Multinomial") comp_srvage_like_vals <- c(comp_srvage_like_vals, 0)
    if(SrvAgeComps_LikeType[f] == "Dirichlet-Multinomial") comp_srvage_like_vals <- c(comp_srvage_like_vals, 1)
    if(SrvAgeComps_LikeType[f] == "iid-Logistic-Normal") comp_srvage_like_vals <- c(comp_srvage_like_vals, 2)
    if(SrvAgeComps_LikeType[f] == "1d-Logistic-Normal") comp_srvage_like_vals <- c(comp_srvage_like_vals, 3)
    if(SrvAgeComps_LikeType[f] == "2d-Logistic-Normal") comp_srvage_like_vals <- c(comp_srvage_like_vals, 4)
    if(SrvAgeComps_LikeType[f] == "3d-Logistic-Normal") comp_srvage_like_vals <- c(comp_srvage_like_vals, 5)
    collect_message(paste("Survey Age Composition Likelihoods", "for survey fleet", f, "specified as:" , SrvAgeComps_LikeType[f]))
  } # end f loop

  # Specifying survey len composition values
  comp_srvlen_like_vals <- vector()
  for(f in 1:input_list$data$n_srv_fleets) {
    if(SrvLenComps_LikeType[f] == 'none') comp_srvlen_like_vals <- c(comp_srvlen_like_vals, 999)
    if(SrvLenComps_LikeType[f] == "Multinomial") comp_srvlen_like_vals <- c(comp_srvlen_like_vals, 0)
    if(SrvLenComps_LikeType[f] == "Dirichlet-Multinomial") comp_srvlen_like_vals <- c(comp_srvlen_like_vals, 1)
    if(SrvLenComps_LikeType[f] == "iid-Logistic-Normal") comp_srvlen_like_vals <- c(comp_srvlen_like_vals, 2)
    if(SrvLenComps_LikeType[f] == "1d-Logistic-Normal") comp_srvlen_like_vals <- c(comp_srvlen_like_vals, 3)
    if(SrvLenComps_LikeType[f] == "2d-Logistic-Normal") comp_srvlen_like_vals <- c(comp_srvlen_like_vals, 4)
    if(SrvLenComps_LikeType[f] == "3d-Logistic-Normal") comp_srvlen_like_vals <- c(comp_srvlen_like_vals, 5)
    collect_message(paste("Survey Length Composition Likelihoods", "for survey fleet", f, "specified as:" , SrvLenComps_LikeType[f]))
  } # end f loop

  # Specifying composition type (ages)
  SrvAgeComps_Type_Mat <- array(NA, dim = c(length(input_list$data$years), input_list$data$n_srv_fleets))
  for(i in 1:length(SrvAgeComps_Type)) {

    # Extract out components from list
    tmp <- SrvAgeComps_Type[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    comps_type_tmp <- tmp_vec[1] # get composition type
    year_range <- as.numeric(unlist(strsplit(tmp_vec[3], "-"))) # get year range
    years <- year_range[1]:year_range[2] # get sequence of years
    fleet <- as.numeric(tmp_vec[5]) # extract fleet index

    # Checking character string
    if(!comps_type_tmp %in% c("agg", "spltRspltS", "spltRjntS", "jntRjntS", 'none')) stop("SrvAgeComps_Type not specified correctly. Must be one of: agg, spltRspltS, spltRjntS, jntRjntS, none")
    if(!fleet %in% c(1:input_list$data$n_srv_fleets)) stop("Invalid fleet specified for SrvAgeComps_Type This needs to be specified as CompType_Year_x-y_Fleet_x")

    # define composition types
    if(comps_type_tmp == "agg") comps_type_val <- 0
    if(comps_type_tmp == "spltRspltS") comps_type_val <- 1
    if(comps_type_tmp == "spltRjntS") comps_type_val <- 2
    if(comps_type_tmp == "jntRjntS") comps_type_val <- 3
    if(comps_type_tmp == "none") comps_type_val <- 999

    # input into matrix
    SrvAgeComps_Type_Mat[years,fleet] <- comps_type_val
  } # end i

  if(any(is.na(SrvAgeComps_Type_Mat))) stop("SrvAgeComps_Type_Mat is returning an NA. Did you update the year range of SrvAgeComps_Type_Mat?")

  SrvLenComps_Type_Mat <- array(NA, dim = c(length(input_list$data$years), input_list$data$n_srv_fleets))
  for(i in 1:length(SrvLenComps_Type)) {

    # Extract out components from list
    tmp <- SrvLenComps_Type[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    comps_type_tmp <- tmp_vec[1] # get composition type
    year_range <- as.numeric(unlist(strsplit(tmp_vec[3], "-"))) # get year range
    years <- year_range[1]:year_range[2] # get sequence of years
    fleet <- as.numeric(tmp_vec[5]) # extract fleet index

    # Checking character string
    if(!comps_type_tmp %in% c("agg", "spltRspltS", "spltRjntS", "jntRjntS", 'none')) stop("SrvLenComps_Type not specified correctly. Must be one of: agg, spltRspltS, spltRjntS, jntRjntS, none")
    if(!fleet %in% c(1:input_list$data$n_srv_fleets)) stop("Invalid fleet specified for SrvLenComps_Type This needs to be specified as CompType_Year_x-y_Fleet_x")

    # define composition types
    if(comps_type_tmp == "agg") comps_type_val <- 0
    if(comps_type_tmp == "spltRspltS") comps_type_val <- 1
    if(comps_type_tmp == "spltRjntS") comps_type_val <- 2
    if(comps_type_tmp == "jntRjntS") comps_type_val <- 3
    if(comps_type_tmp == "none") comps_type_val <- 999

    # input into matrix
    SrvLenComps_Type_Mat[years,fleet] <- comps_type_val
  } # end i

  if(any(is.na(SrvLenComps_Type_Mat))) stop("SrvLenComps_Type_Mat is returning an NA. Did you update the year range of SrvLenComps_Type_Mat?")

  # Input data list
  input_list$data$ObsSrvIdx <- ObsSrvIdx
  input_list$data$ObsSrvIdx_SE <- ObsSrvIdx_SE
  input_list$data$UseSrvIdx <- UseSrvIdx
  input_list$data$ObsSrvAgeComps <- ObsSrvAgeComps
  input_list$data$UseSrvAgeComps <- UseSrvAgeComps
  input_list$data$ObsSrvLenComps <- ObsSrvLenComps
  input_list$data$UseSrvLenComps <- UseSrvLenComps
  input_list$data$SrvAgeComps_LikeType <- comp_srvage_like_vals
  input_list$data$SrvLenComps_LikeType <- comp_srvlen_like_vals
  input_list$data$SrvAgeComps_Type <- SrvAgeComps_Type_Mat
  input_list$data$SrvLenComps_Type <- SrvLenComps_Type_Mat
  input_list$data$srv_idx_type <- srv_idx_type_vals

  # initialize how to aggregate survey age comps
  if(is.null(SrvAge_comp_agg_type)) {
    input_list$data$SrvAge_comp_agg_type <- rep(NA, input_list$data$n_srv_fleets)
  } else input_list$data$SrvAge_comp_agg_type <- SrvAge_comp_agg_type

  # Initialize how to aggregate survey length comps
  if(is.null(SrvLen_comp_agg_type)) {
    input_list$data$SrvLen_comp_agg_type <- rep(NA, input_list$data$n_srv_fleets)
  } else input_list$data$SrvLen_comp_agg_type <- SrvLen_comp_agg_type

  # Setup ISS stuff given differences in how composition data are structured
  if(is.null(ISS_SrvAgeComps)) {
    collect_message("No ISS is specified for SrvAgeComps. ISS weighting is calculated by summing up values from ObsSrvAgeComps each year")
    ISS_SrvAgeComps <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_sexes, input_list$data$n_srv_fleets))
    for(y in 1:length(input_list$data$years)) {
      for(f in 1:input_list$data$n_srv_fleets) {
        # Survey Age Compositions
        # if aggregated across sexes and regions (0) or joint across sexes and regions (3)
        if(input_list$data$SrvAgeComps_Type[y,f] %in% c(0, 3)) ISS_SrvAgeComps[1,y,1,f] <- sum(input_list$data$ObsSrvAgeComps[,y,,,f])
        # if split by region and sex
        if(input_list$data$SrvAgeComps_Type[y,f] == 1) ISS_SrvAgeComps[,y,,f] <- apply(input_list$data$ObsSrvAgeComps[,y,,,f, drop = FALSE], c(1,4), sum)
        # if split by region, joint by sex
        if(input_list$data$SrvAgeComps_Type[y,f] == 2) ISS_SrvAgeComps[,y,1,f] <- apply(input_list$data$ObsSrvAgeComps[,y,,,f, drop = FALSE], 1, sum)
      } # end f loop
    } # end y loop
  }

  if(is.null(ISS_SrvLenComps)) {
    collect_message("No ISS is specified for SrvLenComps. ISS weighting is calculated by summing up values from ObsSrvLenComps each year")
    ISS_SrvLenComps <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_sexes, input_list$data$n_srv_fleets))
    for(y in 1:length(input_list$data$years)) {
      for(f in 1:input_list$data$n_srv_fleets) {
        # Survey Length Compositions
        if(is.null(ISS_SrvLenComps)) {
          # if aggregated across sexes and regions (0) or joint across sexes and regions (3)
          if(input_list$data$SrvLenComps_Type[y,f] %in% c(0, 3)) ISS_SrvLenComps[1,y,1,f] <- sum(input_list$data$ObsSrvLenComps[,y,,,f])
          # if split by region and sex
          if(input_list$data$SrvLenComps_Type[y,f] == 1) ISS_SrvLenComps[,y,,f] <- apply(input_list$data$ObsSrvLenComps[,y,,,f, drop = FALSE], c(1,4), sum)
          # if split by region, joint by sex
          if(input_list$data$SrvLenComps_Type[y,f] == 2) ISS_SrvLenComps[,y,1,f] <- apply(input_list$data$ObsSrvLenComps[,y,,,f, drop = FALSE], 1, sum)
        }
      } # end f loop
    } # end y loop
  }

  # Input into list
  input_list$data$ISS_SrvAgeComps <- ISS_SrvAgeComps
  input_list$data$ISS_SrvLenComps <- ISS_SrvLenComps

  # input parameter list
  starting_values <- list(...)

  # Dispersion parameters for the survey age comps
  if("ln_SrvAge_theta" %in% names(starting_values)) input_list$par$ln_SrvAge_theta <- starting_values$ln_SrvAge_theta
  else input_list$par$ln_SrvAge_theta <- array(0, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_srv_fleets))

  # logistic normal correlation parameters for survey age comps
  if("SrvAge_corr_pars" %in% names(starting_values)) input_list$par$SrvAge_corr_pars <- starting_values$SrvAge_corr_pars
  else input_list$par$SrvAge_corr_pars <- array(0.01, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_srv_fleets, 3))

  # aggregated
  if("ln_SrvAge_theta_agg" %in% names(starting_values)) input_list$par$ln_SrvAge_theta_agg <- starting_values$ln_SrvAge_theta_agg
  else input_list$par$ln_SrvAge_theta_agg <- array(0, dim = c(input_list$data$n_srv_fleets))

  # aggregated correlation parameters
  if("SrvAge_corr_pars_agg" %in% names(starting_values)) input_list$par$SrvAge_corr_pars_agg <- starting_values$SrvAge_corr_pars_agg
  else input_list$par$SrvAge_corr_pars_agg <- array(0.01, dim = c(input_list$data$n_srv_fleets))

  # Dispersion parameters for survey length comps
  if("ln_SrvLen_theta" %in% names(starting_values)) input_list$par$ln_SrvLen_theta <- starting_values$ln_SrvLen_theta
  else input_list$par$ln_SrvLen_theta <- array(0, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_srv_fleets))

  # logistic normal correlation parameters for survey length comps
  if("SrvLen_corr_pars" %in% names(starting_values)) input_list$par$SrvLen_corr_pars <- starting_values$SrvLen_corr_pars
  else input_list$par$SrvLen_corr_pars <- array(0.01, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_srv_fleets, 3))

  # aggregated
  if("ln_SrvLen_theta_agg" %in% names(starting_values)) input_list$par$ln_SrvLen_theta_agg <- starting_values$ln_SrvLen_theta_agg
  else input_list$par$ln_SrvLen_theta_agg <- array(0, dim = c(input_list$data$n_srv_fleets))

  if("SrvLen_corr_pars_agg" %in% names(starting_values)) input_list$par$SrvLen_corr_pars_agg <- starting_values$SrvLen_corr_pars_agg
  else input_list$par$SrvLen_corr_pars_agg <- array(0.01, dim = c(input_list$data$n_srv_fleets))

  # Setup counters
  counter_srvage_agg <- 1
  counter_srvage <- 1
  counter_srvlen_agg <- 1
  counter_srvlen <- 1
  counter_srvage_corr <- 1
  counter_srvage_corr_agg <- 1
  counter_srvlen_corr <- 1
  counter_srvlen_corr_agg <- 1

  # Setup mapping list
  map_SrvAge_theta <- input_list$par$ln_SrvAge_theta # initialize array to set up mapping
  map_SrvLen_theta <- input_list$par$ln_SrvLen_theta # initialize array to set up mapping
  map_SrvAge_theta_agg <- input_list$par$ln_SrvAge_theta_agg # initialize array to set up mapping
  map_SrvLen_theta_agg <- input_list$par$ln_SrvLen_theta_agg # initialize array to set up mapping
  map_SrvAge_corr_pars <- input_list$par$SrvAge_corr_pars # initialize array to set up mapping
  map_SrvAge_corr_pars_agg <- input_list$par$SrvAge_corr_pars_agg # initialize array to set up mapping
  map_SrvLen_corr_pars <- input_list$par$SrvLen_corr_pars # initialize array to set up mapping
  map_SrvLen_corr_pars_agg <- input_list$par$SrvLen_corr_pars_agg # initialize array to set up mapping

  # Initialize these as NAs
  map_SrvAge_theta[] <- NA; map_SrvLen_theta[] <- NA; map_SrvAge_theta_agg[] <- NA; map_SrvLen_theta_agg[] <- NA; map_SrvAge_corr_pars[] <- NA; map_SrvAge_corr_pars_agg[] <- NA; map_SrvLen_corr_pars[] <- NA; map_SrvLen_corr_pars_agg[] <- NA

  for(f in 1:input_list$data$n_srv_fleets) {
    # If we are using a multinomial or there aren't any age comps for a given fleet
    if(input_list$data$SrvAgeComps_LikeType[f] == 0 || sum(input_list$data$UseSrvAgeComps[,,f]) == 0) {
      map_SrvAge_theta[,,f] <- NA
      map_SrvAge_theta_agg[f] <- NA
      map_SrvAge_corr_pars[,,f,] <- NA
      map_SrvAge_corr_pars_agg[f] <- NA
    }

    # If we are using a multinomial or there aren't any lenght comps for a given fleet
    if(input_list$data$SrvLenComps_LikeType[f] == 0 || sum(input_list$data$UseSrvLenComps[,,f]) == 0) {
      map_SrvLen_theta[,,f] <- NA
      map_SrvLen_theta_agg[f] <- NA
      map_SrvLen_corr_pars[,,f,] <- NA
      map_SrvLen_corr_pars_agg[f] <- NA
    }

    # get unique survey age comp types
    srvage_comp_type <- unique(input_list$data$SrvAgeComps_Type[,f])
    srvlen_comp_type <- unique(input_list$data$SrvLenComps_Type[,f])

    # If aggregated (ages)
    if(any(srvage_comp_type == 0) && input_list$data$SrvAgeComps_LikeType[f] != 0) {
      map_SrvAge_theta_agg[f] <- counter_srvage_agg
      counter_srvage_agg <- counter_srvage_agg + 1 # aggregated

      # correlation parameters if any aggregated srvery ages
      if(input_list$data$SrvAgeComps_LikeType[f] == 3) {
        map_SrvAge_corr_pars_agg[f] <- counter_srvage_corr_agg
        counter_srvage_corr_agg <- counter_srvage_corr_agg + 1 # aggregated
      }
    }

    # If aggregated (lengths)
    if(any(srvlen_comp_type == 0) && input_list$data$SrvLenComps_LikeType[f] != 0) {
      map_SrvLen_theta_agg[f] <- counter_srvlen_agg
      counter_srvlen_agg <- counter_srvlen_agg + 1 # aggregated

      # correlation parameters if any aggregated survey lengths
      if(input_list$data$SrvLenComps_LikeType[f] == 3) {
        map_SrvLen_corr_pars_agg[f] <- counter_srvlen_corr_agg
        counter_srvlen_corr_agg <- counter_srvlen_corr_agg + 1 # aggregated
      }
    }

    # joint by sex and region (ages)
    if(any(srvage_comp_type == 3) && input_list$data$SrvAgeComps_LikeType[f] != 0) {
      map_SrvAge_theta[1,1,f] <- counter_srvage
      counter_srvage <- counter_srvage + 1 # joint by sex and region

      # if this is logistic normal, with 3d correlation
      if(input_list$data$SrvAgeComps_LikeType[f] == 5) {
        for(i in 1:3) {
          if(i == 2 && input_list$data$n_sexes == 1) next # skip if we only have 1 sex
          map_SrvAge_corr_pars[1,1,f,i] <- counter_srvage_corr
          counter_srvage_corr <- counter_srvage_corr + 1
        } # end i
      } # end if 3d
    }

    # joint by sex and region (lengths)
    if(any(srvlen_comp_type == 3) && input_list$data$SrvLenComps_LikeType[f] != 0) {
      map_SrvLen_theta[1,1,f] <- counter_srvlen
      counter_srvlen <- counter_srvlen + 1 # joint by sex and region

      # if this is logistic normal, with 3d correlation
      if(input_list$data$SrvLenComps_LikeType[f] == 5) {
        for(i in 1:3) {
          if(i == 2 && input_list$data$n_sexes == 1) next # skip if we only have 1 sex
          map_SrvLen_corr_pars[1,1,f,i] <- counter_srvlen_corr
          counter_srvlen_corr <- counter_srvlen_corr + 1
        } # end i
      } # end if 3d
    }

    # Loop through to make sure mapping stuff off correctly
    for(r in 1:input_list$data$n_regions) {
      for(s in 1:input_list$data$n_sexes) {

        # if split by sex and region (ages)
        if(any(srvage_comp_type == 1) && input_list$data$SrvAgeComps_LikeType[f] != 0) {
          map_SrvAge_theta[r,s,f] <- counter_srvage
          counter_srvage <- counter_srvage + 1 # split by sex and region

          # if logistic normal 1dar1 by age only (unique age correlations for age each region and sex)
          if(input_list$data$SrvAgeComps_LikeType[f] == 3) {
            map_SrvAge_corr_pars[r,s,f,1] <- counter_srvage_corr
            counter_srvage_corr <- counter_srvage_corr + 1
          }
        }

        # if split by sex and region (lengths)
        if(any(srvlen_comp_type == 1) && input_list$data$SrvLenComps_LikeType[f] != 0) {
          map_SrvLen_theta[r,s,f] <- counter_srvlen
          counter_srvlen <- counter_srvlen + 1 # split by sex and region

          # if logistic normal 1dar1 by length only (unique length correlations for length each region and sex)
          if(input_list$data$SrvLenComps_LikeType[f] == 3) {
            map_SrvLen_corr_pars[r,s,f,1] <- counter_srvlen_corr
            counter_srvlen_corr <- counter_srvlen_corr + 1
          }
        }

        # joint by sex, split by region (ages)
        if(any(srvage_comp_type == 2) && input_list$data$SrvAgeComps_LikeType[f] != 0 && s == 1) {
          map_SrvAge_theta[r,1,f] <- counter_srvage
          counter_srvage <- counter_srvage + 1 # joint by sex, split by region

          # if logistic normal, 2d, where 1dar1 by age, constant correlation by sex
          if(input_list$data$SrvAgeComps_LikeType[f] == 4) {
            for(i in 1:2) {
              if(i == 2 && input_list$data$n_sexes == 1) next # skip if we only have 1 sex
              map_SrvAge_corr_pars[r,1,f,i] <- counter_srvage_corr
              counter_srvage_corr <- counter_srvage_corr + 1
            } # end i
          } # end if
        }

        # joint by sex, split by region (lengths)
        if(any(srvlen_comp_type == 2) && input_list$data$SrvLenComps_LikeType[f] != 0 && s == 1) {
          map_SrvLen_theta[r,1,f] <- counter_srvlen
          counter_srvlen <- counter_srvlen + 1 # joint by sex, split by region

          # if logistic normal, 2d, where 1dar1 by length, constant correlation by sex
          if(input_list$data$SrvLenComps_LikeType[f] == 4) {
            for(i in 1:2) {
              if(i == 2 && input_list$data$n_sexes == 1) next # skip if we only have 1 sex
              map_SrvLen_corr_pars[r,1,f,i] <- counter_srvlen_corr
              counter_srvlen_corr <- counter_srvlen_corr + 1
            } # end i
          } # end if
        }

      } # end s loop
    } # end r loop
  } # end f loop

  # Input into mapping list
  input_list$map$ln_SrvAge_theta <- factor(map_SrvAge_theta)
  input_list$map$ln_SrvLen_theta <- factor(map_SrvLen_theta)
  input_list$map$ln_SrvAge_theta_agg <- factor(map_SrvAge_theta_agg)
  input_list$map$ln_SrvLen_theta_agg <- factor(map_SrvLen_theta_agg)
  input_list$map$SrvAge_corr_pars_agg <- factor(map_SrvAge_corr_pars_agg)
  input_list$map$SrvAge_corr_pars <- factor(map_SrvAge_corr_pars)
  input_list$map$SrvLen_corr_pars_agg <- factor(map_SrvLen_corr_pars_agg)
  input_list$map$SrvLen_corr_pars <- factor(map_SrvLen_corr_pars)

  # Print all messages if verbose is TRUE
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}


#' Setup survey selectivity and catchability specifications
#'
#' @param input_list List containing a data list, parameter list, and map list
#' @param srv_sel_blocks Specification for survey selectivity blocks as unique numbers for a given region and fleet, array dimensioned by n_regions, n_years, n_srv_fleets
#' @param srv_sel_model Specification for survey selectivity model for a given region and fleet, array dimensioned by n_regions, n_years, n_srv_fleets, == 0 a50, k, logistic, == 1 gamma dome shaped, == 3, exponential, == 4 a50, a95 logistic
#' @param srv_q_blocks Specification for survey catchability blocks as unique numbers for a given region and fleet, array dimensioned by n_regions, n_years, n_srv_fleets
#' @param ... Additional arguments specifying starting values for survey selectivity and catchability parameters (ln_srv_fixed_sel_pars, ln_srv_q)
#' @param srvsel_pe_pars_spec Specification for survey selectivity process error parameters. If cont_tv_srv_sel is = 0, then this is all fixed and not estimated. Otherwise, the options are: est_all, which estiamtes all parameters, est_shared_r, which estiamtes parameters shared across regions, est_shared_s, which estiamtes parameters shared across sexes, and est_shared_r_s, which estimates these paraemters shared across regions and sexes
#' @param srv_fixed_sel_pars_spec Specification for survey selectivity fixed effects parameters. Options are est_all, which estiamtes all parameters, est_shared_r, which estiamtes parameters shared across regions, est_shared_s, which estiamtes parameters shared across sexes, and est_shared_r_s, which estimates these paraemters shared across regions and sexes
#' @param srv_q_spec Specification for survey catchability. Options are est_all, which estiamtes all parameters across regions, est_shared_r, which estimates parameters shared across regions.
#' @param srv_sel_devs_spec Specificaiton for selectivity process error dviations. Options are est_all, which estimates all deviations, est_shared_r, which shares them across regions, est_shared_s, which shares them across sexes, est_shared_r_s, which shares them across regions and sexes, and est_shared_a, which shares them across age blocks (need to define a number in semipar_age_block_spec), est_shared_r_a, which shares them across regions and age blocks (need to define semipar_age_block_spec), and est_shared_r_a_s, which shares them across regions, ages, and sexes
#' @param semipar_age_block_spec Number (length 1) that needs to be specified if est_shared_a variants are specified. Number represents the number of ages to share deviations for spaced evenly.
#' @param cont_tv_srv_sel Specificaiton for continuous time-varying selectivity, character vector dimensioned by n_srv_fleets, where the character is time variation type, _, Fleet, fleet number. time variation types include (none, iid, rw, 3dmarg, 3dcond, 2dar1), and so if we were to specify iid for fleet 1, this would be iid_Fleet_1.
#' @param corr_opt_semipar Only used if cont_tv_sel is 3,4,5. Allows users to turn off estimation of certain correlation parameters ot be at 0. Options include corr_zero_y, which turns year correlations to 0, corr_zero_a which turns age correaltions to 0, corr_zero_y_a which turns year age correaltions to 0. These options can be used for cont_tv_sel 3,4,5. Additional options include corr_zero_c, which turns cohort correaltions to 0, corr_zero_y_c, which turns cohort and year correaltions to 0, corr_zero_a_c which turns age and cohort correaltions to 0, as well as corr_zero_y_a_c, which turns all correlations to 0, and effectively collapses to iid. These latter options are only available for cont_tv_sel 3,4.
#' @param Use_srv_q_prior Integer specifying whether to use survey q prior or not (0 dont use) (1 use)
#' @param srv_q_prior Survey q priors in normal space, dimensioned by region, block, survey fleet, and 2 (mean, and sd in the 4 dimension of array)
#'
#' @export Setup_Mod_Srvsel_and_Q
#'
Setup_Mod_Srvsel_and_Q <- function(input_list,
                                   cont_tv_srv_sel,
                                   srv_sel_blocks,
                                   srv_sel_model,
                                   Use_srv_q_prior = 0,
                                   srv_q_prior = NA,
                                   srv_q_blocks,
                                   srvsel_pe_pars_spec = NULL,
                                   srv_fixed_sel_pars_spec,
                                   srv_q_spec = NULL,
                                   srv_sel_devs_spec = NULL,
                                   semipar_age_block_spec = NULL,
                                   corr_opt_semipar = NULL,
                                   ...
                                   ) {

  messages_list <<- character(0) # string to attach to for printing messages

  if(!is.null(srvsel_pe_pars_spec)) if(length(srvsel_pe_pars_spec) != input_list$data$n_srv_fleets) stop("srvsel_pe_pars_spec is not length n_srv_fleets")
  if(!is.null(srv_sel_devs_spec)) if(length(srv_sel_devs_spec) != input_list$data$n_srv_fleets) stop("srv_sel_devs_spec is not length n_srv_fleets")
  if(!is.null(corr_opt_semipar)) if(length(corr_opt_semipar) != input_list$data$n_srv_fleets) stop("corr_opt_semipar is not length n_srv_fleets")
  if(!Use_srv_q_prior %in% c(0,1)) stop("Values for Use_srv_q_prior are not valid. They are == 0 (don't use prior), or == 1 (use prior)")
  collect_message("Survey Catchability priors are: ", ifelse(Use_srv_q_prior == 0, "Not Used", "Used"))

  # define for continuous time-varying selectivity
  cont_tv_srv_sel_mat <- array(NA, dim = c(input_list$data$n_regions, input_list$data$n_srv_fleets))
  cont_tv_map <- data.frame(type = c("none", "iid", "rw", "3dmarg", "3dcond", "2dar1"),
                            num = c(0,1,2,3,4,5)) # set up values we map to

  for(i in 1:length(cont_tv_srv_sel)) {
    # Extract out components from list
    tmp <- cont_tv_srv_sel[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    cont_tv_type <- tmp_vec[1] # get continuous selex type
    fleet <- as.numeric(tmp_vec[3]) # extract fleet index
    if(!fleet %in% c(1:input_list$data$n_srv_fleets)) stop("Invalid fleet specified for cont_tv_srv_sel This needs to be specified as timevarytype_Fleet_x")
    if(!cont_tv_type %in% c(cont_tv_map$type)) stop("cont_tv_srv_sel is not correctly specified. This needs to be one of these: none, iid, rw, 3dmarg, 3dcond, 2dar1 (the timevarytypes) and specified as timevarytype_Fleet_x")
    cont_tv_srv_sel_mat[,fleet] <- cont_tv_map$num[which(cont_tv_map$type == cont_tv_type)]
    collect_message("Continuous survey time-varying selectivity specified as: ", cont_tv_type, " for survey fleet ", fleet)
  }

  # define survey time blocks
  srv_sel_blocks_arr = array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets))
  for(i in 1:length(srv_sel_blocks)) {
    # Extract out components from list
    tmp <- srv_sel_blocks[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    if(!tmp_vec[1] %in% c("none", "Block")) stop("Survey Selectivity Blocks not correctly specified. This should be either none_Fleet_x or Block_x_Year_x-y_Fleet_x")
    # extract out fleets if constant
    if(tmp_vec[1] == "none") {
      fleet <- as.numeric(tmp_vec[3]) # get fleet number
      srv_sel_blocks_arr[,,fleet] <- 1 # input only 1 survey time block
    }
    if(tmp_vec[1] == "Block") {
      block_val <- as.numeric(tmp_vec[2]) # get block value
      year_range <- as.numeric(unlist(strsplit(tmp_vec[4], "-"))) # get year range
      years <- year_range[1]:year_range[2] # get sequence of years
      fleet <- as.numeric(tmp_vec[6]) # get fleet number
      srv_sel_blocks_arr[,years,fleet] <- block_val
    }
  }

  for(f in 1:input_list$data$n_srv_fleets) collect_message(paste("Survey Selectivity Time Blocks for survey", f, "is specified at:", length(unique(srv_sel_blocks_arr[,,f]))))

  # Setup survey selectivity models (functional forms)
  sel_map <- data.frame(sel = c('logist1', "gamma", "exponential", "logist2"),
                        num = c(0,1,2,3)) # set up values we can map to
  srv_sel_model_arr = array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets))
  for(i in 1:length(srv_sel_model)) {
    # Extract out components from list
    tmp <- srv_sel_model[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    sel_type <- tmp_vec[1] # get selectivity type
    fleet <- as.numeric(tmp_vec[3]) # get fleet number
    if(!sel_type %in% c(sel_map$sel)) stop("srv_sel_model is not correctly specified. This needs to be one of these: logist1, gamma, exponential, logist2 (the seltypes) and specified as seltype_Fleet_x")
    if(!fleet %in% c(1:input_list$data$n_srv_fleets)) stop("Invalid fleet specified for srv_sel_model This needs to be specified as seltype_Fleet_x")
    srv_sel_model_arr[,,fleet] <- sel_map$num[which(sel_map$sel == sel_type)]
    collect_message("Survey selectivity functional form specified as:", sel_type, " for survey fleet ", f)
  } # end i loop

  if(any(is.na(srv_sel_model_arr))) stop("Survey Selectivity Blocks are returning an NA. Did you update the year range of srv_sel_blocks?")

  # setup survey catchability blocks
  srv_q_blocks_arr = array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets))
  for(i in 1:length(srv_q_blocks)) {
    # Extract out components from list
    tmp <- srv_q_blocks[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    if(!tmp_vec[1] %in% c("none", "Block")) stop("Survey Catchability Blocks not correctly specified. This should be either none_Fleet_x or Block_x_Year_x-y_Fleet_x")
    # extract out fleets if constant
    if(tmp_vec[1] == "none") {
      fleet <- as.numeric(tmp_vec[3]) # get fleet number
      srv_q_blocks_arr[,,fleet] <- 1 # input only 1 survey catchability time block
    }
    if(tmp_vec[1] == "Block") {
      block_val <- as.numeric(tmp_vec[2]) # get block value
      year_range <- as.numeric(unlist(strsplit(tmp_vec[4], "-"))) # get year range
      years <- year_range[1]:year_range[2] # get sequence of years
      fleet <- as.numeric(tmp_vec[6]) # get fleet number
      srv_q_blocks_arr[,years,fleet] <- block_val # input catchability time block
    }
  }

  if(any(is.na(srv_q_blocks_arr))) stop("Survey Catchability Blocks are returning an NA. Did you update the year range of srv_q_blocks?")

  for(f in 1:input_list$data$n_srv_fleets) collect_message(paste("Survey Catchability Time Blocks for fishery", f, "is specified at:", length(unique(srv_q_blocks_arr[,,f]))))

  # Setup data input list
  input_list$data$cont_tv_srv_sel <- cont_tv_srv_sel_mat
  input_list$data$srv_sel_blocks <- srv_sel_blocks_arr
  input_list$data$srv_sel_model <- srv_sel_model_arr
  input_list$data$srv_q_blocks <- srv_q_blocks_arr
  input_list$data$srv_q_prior <- srv_q_prior
  input_list$data$Use_srv_q_prior <- Use_srv_q_prior

  # Set up parameter inputs
  starting_values <- list(...)

  # Survey selectivity process error parameters
  # Survey selectivity fixed effects
  # Figure out number of selectivity parameters for a given functional form
  unique_srvsel_vals <- unique(as.vector(input_list$data$srv_sel_model))
  sel_pars_vec <- vector() # create empty vector to populate

  for(i in 1:length(unique_srvsel_vals)) {
    if(unique_srvsel_vals[i] %in% c(2)) sel_pars_vec[i] <- 1
    if(unique_srvsel_vals[i] %in% c(0,1,3)) sel_pars_vec[i] <- 2
  } # end i loop

  max_srvsel_blks <- max(apply(input_list$data$srv_sel_blocks, c(1,3), FUN = function(x) length(unique(x)))) # figure out maximum number of survey selectivity blocks for a given reigon and fleet
  max_srvsel_pars <- max(sel_pars_vec) # maximum number of selectivity parameters across all forms
  if("ln_srv_fixed_sel_pars" %in% names(starting_values)) input_list$par$ln_srv_fixed_sel_pars <- starting_values$ln_srv_fixed_sel_pars
  else input_list$par$ln_srv_fixed_sel_pars <- array(0, dim = c(input_list$data$n_regions, max_srvsel_pars, max_srvsel_blks, input_list$data$n_sexes, input_list$data$n_srv_fleets))

  # Survey catchability
  max_srvq_blks <- max(apply(input_list$data$srv_q_blocks, c(1,3), FUN = function(x) length(unique(x)))) # figure out maximum number of survey catchability blocks for a given reigon and fleet
  if("ln_srv_q" %in% names(starting_values)) input_list$par$ln_srv_q <- starting_values$ln_srv_q
  else input_list$par$ln_srv_q <- array(0, dim = c(input_list$data$n_regions, max_srvq_blks, input_list$data$n_srv_fleets))

  # Survey selectivity process error parameters
  if("srvsel_pe_pars" %in% names(starting_values)) input_list$par$srvsel_pe_pars <- starting_values$srvsel_pe_pars
  else input_list$par$srvsel_pe_pars <- array(0, dim = c(input_list$data$n_regions, max(max_srvsel_pars, 4), input_list$data$n_sexes, input_list$data$n_srv_fleets)) # dimensioned 4 as the max number of pars for process errors (e.g., sigmas), and then just map off if not using

  # Survey selectivity deviations
  if("ln_srvsel_devs" %in% names(starting_values)) input_list$par$ln_srvsel_devs <- starting_values$ln_srvsel_devs
  else input_list$par$ln_srvsel_devs <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), length(input_list$data$ages), input_list$data$n_sexes, input_list$data$n_srv_fleets))

  # Setup mapping list
  # Initialize counter and mapping array for fixed effects survey selectivity
  srv_fixed_sel_pars_counter <- 1
  map_srv_fixed_sel_pars <- input_list$par$ln_srv_fixed_sel_pars
  map_srv_fixed_sel_pars[] <- NA

  for(f in 1:input_list$data$n_srv_fleets) {
    for(r in 1:input_list$data$n_regions) {

      if(!srv_fixed_sel_pars_spec[f] %in% c("est_all", "est_shared_r", "est_shared_r_s", "fix"))
        stop("srv_fixed_sel_pars_spec not correctly specfied. Should be one of these: est_all, est_shared_r, est_shared_r_s, fix")

      # Figure out max number of selectivity parameters for a given region and fleet
      if(unique(input_list$data$srv_sel_model[r,,f]) %in% 2) max_sel_pars <- 1 # exponential
      if(unique(input_list$data$srv_sel_model[r,,f]) %in% c(0,1,3)) max_sel_pars <- 2 # logistic a50, k, gamma, and logistic a50, a95

      # Extract number of survey selectivity blocks
      srvsel_blocks_tmp <- unique(as.vector(input_list$data$srv_sel_blocks[r,,f]))

      for(s in 1:input_list$data$n_sexes) {
        for(b in 1:length(srvsel_blocks_tmp)) {
          for(i in 1:max_sel_pars) {

            # Estimate all selectivity fixed effects parameters within the constraints of the defined blocks
            if(srv_fixed_sel_pars_spec[f] == "est_all") {
              map_srv_fixed_sel_pars[r,i,b,s,f] <- srv_fixed_sel_pars_counter
              srv_fixed_sel_pars_counter <- srv_fixed_sel_pars_counter + 1
            } # end if

            # Estimating process error parameters shared across regions (but unique for each sex, fleet, parameter)
            if(srv_fixed_sel_pars_spec[f] == 'est_shared_r' && r == 1) {
              for(rr in 1:input_list$data$n_regions) {
                if(srvsel_blocks_tmp[b] %in% input_list$data$srv_sel_blocks[rr,,f]) {
                  map_srv_fixed_sel_pars[rr, i, b, s, f] <- srv_fixed_sel_pars_counter
                } # end if
              } # end rr loop
              srv_fixed_sel_pars_counter <- srv_fixed_sel_pars_counter + 1
            } # end if

            # Estimating process error parameters shared across sexes (but unique for each region, fleet, parameter)
            if(srv_fixed_sel_pars_spec[f] == 'est_shared_s' && s == 1) {
              for(ss in 1:input_list$data$n_sexes) {
                if(srvsel_blocks_tmp[b] %in% input_list$data$srv_sel_blocks[r,,f]) {
                  map_srv_fixed_sel_pars[r, i, b, ss, f] <- srv_fixed_sel_pars_counter
                } # end if
              } # end ss loop
              srv_fixed_sel_pars_counter <- srv_fixed_sel_pars_counter + 1
            } # end if

            # Estimating process error parameters shared across regions and sexes (but unique for each fleet, parameter)
            if(srv_fixed_sel_pars_spec[f] == 'est_shared_r_s' && r == 1 && s == 1) {
              for(rr in 1:input_list$data$n_regions) {
                for(ss in 1:input_list$data$n_sexes) {
                  if(srvsel_blocks_tmp[b] %in% input_list$data$srv_sel_blocks[rr,,f]) {
                    map_srv_fixed_sel_pars[rr, i, b, ss, f] <- srv_fixed_sel_pars_counter
                  } # end if
                } # end ss loop
              } # end rr loop
              srv_fixed_sel_pars_counter <- srv_fixed_sel_pars_counter + 1
            } # end if

          } # end i loop
        } # end b loop
      } # end s loop
    } # end r loop

    # fix all parameters
    if(srv_fixed_sel_pars_spec[f] == "fix") map_srv_fixed_sel_pars[,,,,f] <- NA
    collect_message("srv_fixed_sel_pars_spec is specified as: ", srv_fixed_sel_pars_spec[f], " for survey fleet ", f)

  } # end f loop

  # Initialize counter and mapping array for survey catchability
  srv_q_counter <- 1
  map_srv_q <- input_list$par$ln_srv_q
  map_srv_q[] <- NA

  for(f in 1:input_list$data$n_srv_fleets) {

    if(!is.null(srv_q_spec)) {
      if(!srv_q_spec[f] %in% c("est_all", "est_shared_r", "fix"))
        stop("srv_q_spec not correctly specfied. Should be one of these: est_all, est_shared_r, fix")
    }

    for(r in 1:input_list$data$n_regions) {
      # Extract number of survey catchability and region blocks
      srvq_blocks_tmp <- unique(as.vector(input_list$data$srv_q_blocks[r,,f]))
      if(sum(input_list$data$UseSrvIdx[r,,f]) == 0) {
        map_srv_q[r,,f] <- NA # fix parameters if we are not using survey indices for these fleets and regions
      } else {
        for(b in 1:length(srvq_blocks_tmp)) {
          # Estimate for all regions
          if(srv_q_spec[f] == 'est_all') {
            map_srv_q[r,b,f] <- srv_q_counter
            srv_q_counter <- srv_q_counter + 1
          }
          # Estimate but share q across regions
          if(srv_q_spec[f] == 'est_shared_r' && r == 1) {
            for(rr in 1:input_list$data$n_regions) {
              if(srvq_blocks_tmp[b] %in% input_list$data$srv_q_blocks[rr,,f]) {
                map_srv_q[rr, b, f] <- srv_q_counter
              } # end if
            } # end rr loop
            srv_q_counter <- srv_q_counter + 1
          } # end if
        } # end b loop

        # fix all parameters
        if(srv_q_spec[f] == 'fix') map_srv_q[,,f] <- NA
      } # end else loop
    } # end r loop
    collect_message("srv_q_spec is specified as: ", srv_q_spec[f], " for survey fleet ", f)
  } # end f loop

  # Initialize counter and mapping array for survey process errors
  srvsel_pe_pars_counter <- 1 # initalize counter
  map_srvsel_pe_pars <- input_list$par$srvsel_pe_pars # initalize array
  map_srvsel_pe_pars[] <- NA

  # survey process error parameters
  for(f in 1:input_list$data$n_srv_fleets) {

    if(!is.null(srvsel_pe_pars_spec)) {
      if(!srvsel_pe_pars_spec[f] %in% c("est_all", "est_shared_r", NA, "est_shared_s", "est_shared_r_s"))
        stop("srvsel_pe_pars_spec not correctly specfied. Should be one of these: est_all, est_shared_r, est_shared_r_s, fix, NA")
    }

    for(r in 1:input_list$data$n_regions) {

      # Figure out max number of selectivity parameters for a given region and fleet
      if(unique(input_list$data$srv_sel_model[r,,f]) %in% 2) max_sel_pars <- 1 # exponential
      if(unique(input_list$data$srv_sel_model[r,,f]) %in% c(0,1,3)) max_sel_pars <- 2 # logistic a50, k, gamma, and logistic a50, a95

      # if no time-variation, then fix all parameters for this fleet
      if(input_list$data$cont_tv_srv_sel[r,f] == 0) {
        map_srvsel_pe_pars[r,,,f] <- NA
      } else { # if we have time-variation
        for(s in 1:input_list$data$n_sexes) {
          # If iid time-variation for this fleet
          if(input_list$data$cont_tv_srv_sel[r,f] == 1) {
            for(i in 1:max_sel_pars) {
              # Estimating all parameters separately (unique for each region, sex, fleet, parameter)
              if(srvsel_pe_pars_spec[f] == "est_all") {
                map_srvsel_pe_pars[r,i,s,f] <- srvsel_pe_pars_counter
                srvsel_pe_pars_counter <- srvsel_pe_pars_counter + 1
              } # end est_all
              # Estimating process error parameters shared across regions (but unique for each sex, fleet, parameter)
              if(srvsel_pe_pars_spec[f] == 'est_shared_r' && r == 1) {
                map_srvsel_pe_pars[,i,s,f] <- srvsel_pe_pars_counter
                srvsel_pe_pars_counter <- srvsel_pe_pars_counter + 1
              }
              # Estimating process error parameters shared across sexes (but unique for each region, fleet, parameter)
              if(srvsel_pe_pars_spec[f] == 'est_shared_s' && s == 1) {
                map_srvsel_pe_pars[r,i,,f] <- srvsel_pe_pars_counter
                srvsel_pe_pars_counter <- srvsel_pe_pars_counter + 1
              }
              # Estimating process error parameters shared across regions and sexes (but unique for each fleet, parameter)
              if(srvsel_pe_pars_spec[f] == 'est_shared_r_s' && r == 1 && s == 1) {
                map_srvsel_pe_pars[,i,,f] <- srvsel_pe_pars_counter
                srvsel_pe_pars_counter <- srvsel_pe_pars_counter + 1
              }
            } # end i loop
          } # end iid variation

          # If random walk time-variation for this fleet
          if(y > 1) {
            if(input_list$data$cont_tv_srv_sel[r,f] == 2) {
              for(i in 1:max_sel_pars) {
                # Estimating all parameters separately (unique for each region, sex, fleet, parameter)
                if(srvsel_pe_pars_spec[f] == "est_all") {
                  map_srvsel_pe_pars[r,i,s,f] <- srvsel_pe_pars_counter
                  srvsel_pe_pars_counter <- srvsel_pe_pars_counter + 1
                } # end est_all
                # Estimating process error parameters shared across regions (but unique for each sex, fleet, parameter)
                if(srvsel_pe_pars_spec[f] == 'est_shared_r' && r == 1) {
                  map_srvsel_pe_pars[,i,s,f] <- srvsel_pe_pars_counter
                  srvsel_pe_pars_counter <- srvsel_pe_pars_counter + 1
                }
                # Estimating process error parameters shared across sexes (but unique for each region, fleet, parameter)
                if(srvsel_pe_pars_spec[f] == 'est_shared_s' && s == 1) {
                  map_srvsel_pe_pars[r,i,,f] <- srvsel_pe_pars_counter
                  srvsel_pe_pars_counter <- srvsel_pe_pars_counter + 1
                }
                # Estimating process error parameters shared across regions and sexes (but unique for each fleet, parameter)
                if(srvsel_pe_pars_spec[f] == 'est_shared_r_s' && r == 1 && s == 1) {
                  map_srvsel_pe_pars[,i,,f] <- srvsel_pe_pars_counter
                  srvsel_pe_pars_counter <- srvsel_pe_pars_counter + 1
                }
              } # end i loop
            } # end random walk variation
          } # only estimate if y > 1, otherwise devs set to zero

          # If 3d gmrf or 2dar1
          if(input_list$data$cont_tv_srv_sel[r,f] %in% c(3,4,5)) {

            # Set up indexing to loop through
            if(input_list$data$cont_tv_srv_sel[r,f] %in% c(3,4)) idx = 1:4 # 3dgmrf (1 = pcorr_age, 2 = pcorr_year, 3= pcorr_cohort, 4 = log_sigma)
            if(input_list$data$cont_tv_srv_sel[r,f] %in% c(5)) idx = c(1,2,4) # 2dar1 (1 = pcorr_age, 2 = pcorr_year, 4 = log_sigma)

            for(i in idx) {
              # Estimating all process error parameters
              if(srvsel_pe_pars_spec[f] == "est_all") {
                map_srvsel_pe_pars[r,i,s,f] <- srvsel_pe_pars_counter
                srvsel_pe_pars_counter <- srvsel_pe_pars_counter + 1
              } # end est_all
              # Estimating process error parameters shared across regions (but unique for each sex, fleet, parameter)
              if(srvsel_pe_pars_spec[f] == 'est_shared_r' && r == 1) {
                map_srvsel_pe_pars[,i,s,f] <- srvsel_pe_pars_counter
                srvsel_pe_pars_counter <- srvsel_pe_pars_counter + 1
              }
              # Estimating process error parameters shared across sexes (but unique for each region, fleet, parameter)
              if(srvsel_pe_pars_spec[f] == 'est_shared_s' && s == 1) {
                map_srvsel_pe_pars[r,i,,f] <- srvsel_pe_pars_counter
                srvsel_pe_pars_counter <- srvsel_pe_pars_counter + 1
              }
              # Estimating process error parameters shared across regions and sexes (but unique for each fleet, parameter)
              if(srvsel_pe_pars_spec[f] == 'est_shared_r_s' && r == 1 && s == 1) {
                map_srvsel_pe_pars[,i,,f] <- srvsel_pe_pars_counter
                srvsel_pe_pars_counter <- srvsel_pe_pars_counter + 1
              }
            } # end i loop

            # Options to set correaltions to 0 for 3dgmrf
            if(!is.null(corr_opt_semipar)) {

              if(!corr_opt_semipar[f] %in% c(NA, "corr_zero_y", "corr_zero_a", "corr_zero_y_a", "corr_zero_c", "corr_zero_y_c", "corr_zero_a_c", "corr_zero_y_a_c"))
                stop("corr_opt_semipar not correctly specfied. Should be one of these: NA, corr_zero_y, corr_zero_a, corr_zero_y_a, corr_zero_c, corr_zero_y_c, corr_zero_a_c, corr_zero_y_a_c")

              # if either 2dar1 or 3dgmrf
              if(input_list$data$cont_tv_srv_sel[r,f] %in% c(3,4,5)) {
                if(corr_opt_semipar[f] == "corr_zero_y") map_srvsel_pe_pars[,2,,1] <- NA
                if(corr_opt_semipar[f] == "corr_zero_a") map_srvsel_pe_pars[,1,,1] <- NA
                if(corr_opt_semipar[f] == "corr_zero_y_a") map_srvsel_pe_pars[,1:2,,1] <- NA
              }

              # if 3dgmrf only
              if(input_list$data$cont_tv_srv_sel[r,f] %in% c(3,4)) {
                if(corr_opt_semipar[f] == "corr_zero_c") map_srvsel_pe_pars[,3,,1] <- NA
                if(corr_opt_semipar[f] == "corr_zero_y_c") map_srvsel_pe_pars[,2:3,,1] <- NA
                if(corr_opt_semipar[f] == "corr_zero_a_c") map_srvsel_pe_pars[,c(1,3),,1] <- NA
                if(corr_opt_semipar[f] == "corr_zero_y_a_c") map_srvsel_pe_pars[,1:3,,1] <- NA
              }

              # Reset numbering for mapping off correlation parameters for clarity
              non_na_positions <- which(!is.na(map_srvsel_pe_pars))
              map_srvsel_pe_pars[non_na_positions] <- seq_along(non_na_positions)
              collect_message("corr_opt_semipar is specified as: ", corr_opt_semipar[f], "for survey fleet", f)

            }

          } # end if 3d gmrf marginal or conditional variance

          # fix all parameters
          if(srvsel_pe_pars_spec[f] == "fix") map_srvsel_pe_pars[r,,s,f] <- NA

        } # end s loop
      } # end else
    } # end r loop

    if(!is.null(srvsel_pe_pars_spec)) collect_message("srvsel_pe_pars_spec is specified as: ", srvsel_pe_pars_spec[f], " for survey fleet ", f)

  } # end f loop

  # Initialize counter and mapping array for survey selectivity deviations
  srvsel_devs_counter <- 1
  map_srvsel_devs <- input_list$par$ln_srvsel_devs
  map_srvsel_devs[] <- NA

  for(r in 1:input_list$data$n_regions) {
    for(f in 1:input_list$data$n_srv_fleets) {

      if(!is.null(srv_sel_devs_spec)) {
        if(!srv_sel_devs_spec[f] %in% c(NA, "est_all", "est_shared_r", "est_shared_s", "est_shared_r_s", "est_shared_a", "est_shared_r_a", "est_shared_a_s", "est_shared_r_a_s"))
          stop("srv_sel_devs_spec not correctly specfied. Should be one of these: est_all, est_shared_r, est_shared_r_s, est_shared_a, est_shared_r_a, est_shared_a_s, est_shared_r_a_s")
      }

      # Figure out max number of selectivity parameters for a given region and fleet
      if(unique(input_list$data$srv_sel_model[r,,f]) %in% 2) max_sel_pars <- 1 # exponential
      if(unique(input_list$data$srv_sel_model[r,,f]) %in% c(0,1,3)) max_sel_pars <- 2 # logistic a50, k, gamma, and logistic a50, a95

      for(s in 1:input_list$data$n_sexes) {
        for(y in 1:length(input_list$data$years)) {

          # if no time-variation, then fix all parameters for this fleet
          if(input_list$data$cont_tv_srv_sel[r,f] == 0) {
            map_srvsel_devs[r,y,,s,f] <- NA
          } else {
            # If iid time-variation or random walk for this fleet
            if(input_list$data$cont_tv_srv_sel[r,f] %in% c(1,2)) {
              for(i in 1:max_sel_pars) {
                # Estimating all selectivity deviations across regions, sexes, fleets, and parameter
                if(srv_sel_devs_spec[f] == 'est_all') {
                  map_srvsel_devs[r,y,i,s,f] <- srvsel_devs_counter
                  srvsel_devs_counter <- srvsel_devs_counter + 1
                }
                # Estimating selectivity deviations across sexes, fleets, and parameters, but shared across regions
                if(srv_sel_devs_spec[f] == 'est_shared_r' && r == 1) {
                  map_srvsel_devs[,y,i,s,f] <- srvsel_devs_counter
                  srvsel_devs_counter <- srvsel_devs_counter + 1
                }
                # Estimating selectivity deviations across regions, fleets, and parameters, but shared across sexes
                if(srv_sel_devs_spec[f] == 'est_shared_s' && r == 1) {
                  map_srvsel_devs[r,y,i,,f] <- srvsel_devs_counter
                  srvsel_devs_counter <- srvsel_devs_counter + 1
                }
                # Estimating selectivity deviations across fleets, and parameters, but shared across sexes and regions
                if(srv_sel_devs_spec[f] == 'est_shared_s' && r == 1 && s == 1) {
                  map_srvsel_devs[,y,i,,f] <- srvsel_devs_counter
                  srvsel_devs_counter <- srvsel_devs_counter + 1
                }
              } # end i loop
            } # end iid or random walk variation

            # If 3d gmrf for this fleet
            if(input_list$data$cont_tv_srv_sel[r,f] %in% c(3,4,5)) {
              for(i in 1:length(input_list$data$ages)) {
                # Estimating all selectivity deviations across regions, years and ages (also cohorts baked in year x age)
                if(srv_sel_devs_spec[f] == 'est_all') {
                  map_srvsel_devs[r,y,i,s,f] <- srvsel_devs_counter
                  srvsel_devs_counter <- srvsel_devs_counter + 1
                }
                # Estimating all selectivity deviations across years and ages (also cohorts baked in year x age), but shared across regions
                if(srv_sel_devs_spec[f] == 'est_shared_r' && r == 1) {
                  map_srvsel_devs[,y,i,s,f] <- srvsel_devs_counter
                  srvsel_devs_counter <- srvsel_devs_counter + 1
                }
              } # end i loop

              # Estimating selectivity deviations across years, but blocking structure / pooling across ages
              if(srv_sel_devs_spec[f] == 'est_shared_a') {
                age_block_index <- ceiling(input_list$data$ages / semipar_age_block_spec) # get age block indexing
                for(i in 1:length(unique(age_block_index))) {
                  map_srvsel_devs[r,y,which(age_block_index == i),s,f] <- srvsel_devs_counter
                  srvsel_devs_counter <- srvsel_devs_counter + 1
                } # end i loop
              } # end if sharing age blocks

              # Estimating selectivity deviations across years, but blocking structure / pooling across ages and sharing across regions
              if(srv_sel_devs_spec[f] == 'est_shared_r_a' && r == 1) {
                age_block_index <- ceiling(input_list$data$ages / semipar_age_block_spec) # get age block indexing
                for(i in 1:length(unique(age_block_index))) {
                  map_srvsel_devs[,y,which(age_block_index == i),s,f] <- srvsel_devs_counter
                  srvsel_devs_counter <- srvsel_devs_counter + 1
                } # end i loop
              } # end if sharing age blocks

              # Estimating selectivity deviations across years, but blocking and sharing across ages, and sharing across sexes
              if(srv_sel_devs_spec[f] == 'est_shared_a_s' && s == 1) {
                age_block_index <- ceiling(input_list$data$ages / semipar_age_block_spec) # get age block indexing
                for(i in 1:length(unique(age_block_index))) {
                  map_srvsel_devs[r,y,which(age_block_index == i),,f] <- srvsel_devs_counter
                  srvsel_devs_counter <- srvsel_devs_counter + 1
                } # end i loop
              } # end if sharing age blocks and sexes

              # Estimating selectivity deviations across years, but blocking and sharing across ages, and sharing across sexes  and regions
              if(srv_sel_devs_spec[f] == 'est_shared_r_a_s' && r == 1 && s == 1) {
                age_block_index <- ceiling(input_list$data$ages / semipar_age_block_spec) # get age block indexing
                for(i in 1:length(unique(age_block_index))) {
                  map_srvsel_devs[,y,which(age_block_index == i),,f] <- srvsel_devs_counter
                  srvsel_devs_counter <- srvsel_devs_counter + 1
                } # end i loop
              } # end if sharing age blocks and sexes and regions
            } # end 3d gmrf

          } # end else
        } # end y loop
      } # end s loop

      if(!is.null(srv_sel_devs_spec)) collect_message("srv_sel_devs_spec is specified as: ", srv_sel_devs_spec[f], " for survey fleet ", f, " and region ", r)

    } # end f loop
  } # end r loop

  # Input into list
  input_list$map$srvsel_pe_pars <- factor(map_srvsel_pe_pars)
  input_list$map$ln_srv_fixed_sel_pars <- factor(map_srv_fixed_sel_pars)
  input_list$map$ln_srv_q <- factor(map_srv_q)
  input_list$map$ln_srvsel_devs <- factor(map_srvsel_devs)
  input_list$data$map_ln_srvsel_devs <- array(as.numeric(input_list$map$ln_srvsel_devs), dim = dim(input_list$par$ln_srvsel_devs))
  input_list$data$map_srv_q <- map_srv_q

  # Checking whether survey q dimensions are correct
  if(!is.na(srv_q_prior) || Use_srv_q_prior == 1) if(sum(dim(srv_q_prior) == c(dim(map_srv_q), 2)) != 4) stop("Survey catchability dimensions are not correct. Should be n_regions, max n_blocks, n_srv_fleets, and 2 (where 2 represents the 2 prior parameters - the mean and sd). You can input an NA if not availiable for certain regions or fleets")

  # Print all messages if verbose is TRUE
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}
