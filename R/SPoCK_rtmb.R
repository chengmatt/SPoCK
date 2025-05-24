# version 1 - 12/22/24 (M.LH Cheng)
# Bridge model 23.5 from ADMB to RTMB
# Changed code to be more modular, accommodating any number of fishery and survey fleets
# Rectified errors in fitting to length composition data (normalize proportions at length after conversion from age-length matrix)
# Changed survey composition data to be calculated using survival midyear
# Added options for continuous time-varying selectivity
# Added options for TMB / R-like likelihoods (e.g., dnorm) instead of custom likelihoods
# Added options for dirichlet multinomial likelihood

# version 2 - 12/23/24 (M.LH Cheng)
# Incorporated options to fit age and length composition data as sex-aggregated, split by sex (no sex ratio
# information), and jointly by sex (implicit sex ratio information)
# Added in option for Dirichlet Multinomial likelihood

# version 3 - (M.LH Cheng)
# Coded in spatial dimensions
# Parameters (mean recruitment, recruitment devs, initial age devs,
# selectivity, composition likelihood parameters, catchability,
# mean fishing mortality, fishing mortality deviates) can be estimated spatially
# Incorporated options to allow for estimation of movement parameters across years, ages, and sexes
# Tag integrated model incorporated using a Brownie Tag Attrition Model
# Tag Reporting Rates, Tag Shedding, and Tag Induced Mortality are parameters that can be estimated
# Beta priors for tag reporting rates, dirichlet priors for movement rates
# Incorporated iid, random walk, 2d and 3d correaltions for fishery and survey selectivity
# Added in options for Logistic Normal likelihood

#' Generalized RTMB model
#'
#' @param pars Parameter List
#' @param data Data List
#' @import RTMB
#' @keywords internal
SPoCK_rtmb = function(pars, data) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  RTMB::getAll(pars, data) # load in starting values and data

  # Model Set Up (Containers) -----------------------------------------------
  n_ages = length(ages) # number of ages
  n_yrs = length(years) # number of years
  n_lens = length(lens) # number of lengths

  # Recruitment stuff
  n_est_rec_devs = dim(ln_RecDevs)[2] # number of recruitment deviates estimated
  Rec = array(0, dim = c(n_regions, n_yrs)) # Recruitment
  R0 = rep(0, n_regions) # R0 or mean recruitment

  # Population Dynamics
  init_iter = n_ages * 5 # Number of times to iterate to equilibrium when movement occurs
  Init_NAA = array(0, dim = c(n_regions, n_ages, n_sexes)) # initial age structure
  Init_NAA_next_year = Init_NAA # initial age structure
  NAA = array(data = 0, dim = c(n_regions, n_yrs + 1, n_ages, n_sexes)) # Numbers at age
  ZAA = array(data = 0, dim = c(n_regions, n_yrs, n_ages, n_sexes)) # Total mortality at age
  SAA_mid = array(data = 0, dim = c(n_regions, n_yrs, n_ages, n_sexes)) # Survival at age (midpoint of the year)
  natmort = array(data = 0, dim = c(n_regions, n_yrs, n_ages, n_sexes)) # natural mortality at age
  Total_Biom = array(0, dim = c(n_regions, n_yrs)) # Total biomass
  SSB = array(0, dim = c(n_regions, n_yrs)) # Spawning stock biomass

  # Movement Stuff
  Movement = array(data = 0, dim = c(n_regions, n_regions, n_yrs, n_ages, n_sexes)) # movement "matrix"
  n_move_age_tag_pool = length(move_age_tag_pool) # number of ages to pool for tagging data
  n_move_sex_tag_pool = length(move_sex_tag_pool) # number of sexes to pool for tagging data

  # Tagging Stuff
  Tags_Avail = array(data = 0, dim = c(max_tag_liberty + 1, n_tag_cohorts, n_regions, n_ages, n_sexes)) # Tags availiable for recapture
  Tag_Reporting = array(data = 0, dim = c(n_regions, n_yrs)) # Tag reporting rate
  Pred_Tag_Recap = array(data = 0, dim = c(max_tag_liberty, n_tag_cohorts, n_regions, n_ages, n_sexes)) # predicted recaptures

  # Fishery Processes
  init_F = init_F_prop * exp(ln_F_mean[1]) # initial F for age structure
  Fmort = array(0, dim = c(n_regions, n_yrs, n_fish_fleets)) # Fishing mortality scalar
  FAA = array(data = 0, dim = c(n_regions, n_yrs, n_ages, n_sexes, n_fish_fleets)) # Fishing mortality at age
  CAA = array(data = 0, dim = c(n_regions, n_yrs, n_ages, n_sexes, n_fish_fleets)) # Catch at age
  CAL = array(data = 0, dim = c(n_regions, n_yrs, n_lens, n_sexes, n_fish_fleets)) # Catch at length
  PredCatch = array(0, dim = c(n_regions, n_yrs, n_fish_fleets)) # Predicted catch in weight
  PredFishIdx = array(0, dim = c(n_regions, n_yrs, n_fish_fleets)) # Predicted fishery index
  fish_sel = array(data = 0, dim = c(n_regions, n_yrs, n_ages, n_sexes, n_fish_fleets)) # Fishery selectivity
  fish_q = array(0, dim = c(n_regions, n_yrs, n_fish_fleets)) # Fishery catchability

  # Survey Processes
  SrvIAA = array(data = 0, dim = c(n_regions, n_yrs, n_ages, n_sexes, n_srv_fleets)) # Survey index at age
  SrvIAL = array(data = 0, dim = c(n_regions, n_yrs, n_lens, n_sexes, n_srv_fleets)) # Survey index at length
  PredSrvIdx = array(0, dim = c(n_regions, n_yrs, n_srv_fleets)) # Predicted survey index
  srv_sel = array(data = 0, dim = c(n_regions, n_yrs, n_ages, n_sexes, n_srv_fleets)) # Survey selectivity
  srv_q = array(0, dim = c(n_regions, n_yrs, n_srv_fleets)) # Survey catchability

  # Likelihoods
  Catch_nLL = array(0, dim = c(n_regions, n_yrs, n_fish_fleets)) # Fishery Catch Likelihoods
  FishIdx_nLL = array(0, dim = c(n_regions, n_yrs, n_fish_fleets)) # Fishery Index Likelihoods
  FishAgeComps_nLL = array(data = 0, dim = c(n_regions, n_yrs, n_sexes, n_fish_fleets)) # Fishery Age Comps Likelihoods
  FishLenComps_nLL = array(data = 0, dim = c(n_regions, n_yrs, n_sexes, n_fish_fleets)) # Fishery Length Comps Likelihoods
  SrvIdx_nLL = array(0, dim = c(n_regions, n_yrs, n_srv_fleets)) # Survey Index Likelihoods
  SrvAgeComps_nLL = array(data = 0, dim = c(n_regions, n_yrs, n_sexes, n_srv_fleets)) # Survey Age Comps Likelihoods
  SrvLenComps_nLL = array(data = 0, dim = c(n_regions, n_yrs, n_sexes, n_srv_fleets)) # Survey Length Comps Likelihoods
  Tag_nLL = array(data = 0, dim = c(max_tag_liberty, n_tag_cohorts, n_regions, n_ages, n_sexes)) # Tagging Likelihoods

  # Penalties and Priors
  Fmort_nLL = array(0, dim = c(n_regions, n_yrs, n_fish_fleets)) # Fishing Mortality Deviation penalty
  Rec_nLL = array(0, dim = c(n_regions, n_yrs)) # Recruitment penalty
  Init_Rec_nLL = array(0, dim = c(n_regions, n_ages - 2)) # Initial Recruitment penalty
  bias_ramp = rep(0, n_yrs) # bias ramp from Methot and Taylor 2011
  sel_nLL = 0 # Penalty for selectivity deviations
  fish_q_nLL = 0 # Prior/penalty for fishery q
  srv_q_nLL = 0 # Prior/penalty for survey q
  M_nLL = 0 # Penalty/Prior for natural mortality
  h_nLL = 0 # Prior for steepness
  Movement_nLL = 0 # Penalty for movement rates
  TagRep_nLL = 0 # penalty for tag reporting rate
  jnLL = 0 # Joint negative log likelihood

  # Model Process Equations -------------------------------------------------
  ## Movement Parameters (Set up) --------------------------------------------
  ref_region = 1 # Set up reference region (always set at 0)
  for(r in 1:n_regions) {
    for(y in 1:n_yrs) {
      for(a in 1:n_ages) {
        for(s in 1:n_sexes) {
          move_tmp = rep(0, n_regions) # temporary movement vector to store values (from - to)
          counter = 1  # counter
          for(rr in 1:n_regions) {
            if(rr != ref_region) {
              move_tmp[rr] = move_pars[r,counter,y,a,s] + logit_move_devs[r,counter,y,a]
              counter = counter + 1
            } # end if not reference region
          } # end rr loop
          if(use_fixed_movement == 0) Movement[r,,y,a,s] = exp(move_tmp) / sum(exp(move_tmp)) # multinomial logit transform (basically a softmax) - estimated movement
          if(use_fixed_movement == 1) Movement[r,,y,a,s] = Fixed_Movement[r,,y,a,s] # fixed movement matrix
        } # end s loop
      } # end a loop
    } # end y loop
  } # end r loop

  ## Fishery Selectivity -----------------------------------------------------
  for(r in 1:n_regions) {
    for(y in 1:n_yrs) {
      for(f in 1:n_fish_fleets) {
        fish_sel_blk_idx = fish_sel_blocks[r,y,f] # Get fishery selectivity block index
        for(s in 1:n_sexes) {

          # Extract out fixed-effect selectivity parameters
          tmp_fish_sel_vec = ln_fish_fixed_sel_pars[r,,fish_sel_blk_idx,s,f]

          # Calculate selectivity
          fish_sel[r,y,,s,f] = Get_Selex(Selex_Model = fish_sel_model[r,y,f], # selectivity model
                                         TimeVary_Model = cont_tv_fish_sel[r,f], # time varying model
                                         ln_Pars = tmp_fish_sel_vec, # fixed effect selectivity parameters
                                         ln_seldevs = ln_fishsel_devs[,,,,f, drop = FALSE], # list object to incorporate different dimensions of deviations
                                         Region = r, # region index
                                         Year = y, # year index
                                         Age = ages, # age vector
                                         Sex = s) # sex index

        } # end s loop
      } # end f loop
    } # end y loop
  } # end r loop


  ## Survey Selectivity ------------------------------------------------------
  for(r in 1:n_regions) {
    for(y in 1:n_yrs) {
      for(sf in 1:n_srv_fleets) {
        srv_sel_blk_idx = srv_sel_blocks[r,y,sf] # Get survey selectivity block index
        for(s in 1:n_sexes) {

          # extract temporary selectivity parameters
          tmp_srv_sel_vec = ln_srv_fixed_sel_pars[r,,srv_sel_blk_idx,s,sf]

          # Calculate selectivity
          srv_sel[r,y,,s,sf] = Get_Selex(Selex_Model = srv_sel_model[r,y,sf], # selectivity model
                                         TimeVary_Model = cont_tv_srv_sel[r,sf], # time varying model
                                         ln_seldevs = ln_srvsel_devs[,,,,sf, drop = FALSE], # deviations
                                         ln_Pars = tmp_srv_sel_vec,
                                         Region = r, # region index
                                         Year = y, # year index
                                         Age = ages, # age vector
                                         Sex = s) # sex index

        } # end s loop
      } # end sf loop
    } # end y loop
  } # end r loop


  ## Mortality ---------------------------------------------------------------
  for(r in 1:n_regions) {
    for(y in 1:n_yrs) {

      # Fishing Mortality at Age calculations
      for(f in 1:n_fish_fleets) {
        if(UseCatch[r,y,f] == 0) {
          Fmort[r,y,f] = 0 # Set F to zero when no catch data
        } else {
          if(Catch_Type[y,f] == 0 && est_all_regional_F == 0) {
            Fmort[r,y,f] = exp(ln_F_mean_AggCatch[f] + ln_F_devs_AggCatch[y,f]) # If catch is aggregated across regions
          } else Fmort[r,y,f] = exp(ln_F_mean[r,f] + ln_F_devs[r,y,f]) # Fully selected F
        }
        FAA[r,y,,,f] = Fmort[r,y,f] * fish_sel[r,y,,,f,drop = FALSE] # Fishing mortality at age
      } # f loop

      # Population Mortality and Survival
      natmort[r,y,,1] = exp(ln_M) # get natural mortality (females or single-sex)
      if(n_sexes == 2) natmort[r,y,,2] = exp(ln_M) + M_offset # natural mortality with offset (males)
      ZAA[r,y,,] = apply(FAA[r,y,,,,drop = FALSE],3:4,sum) + natmort[r,y,,] # Total Mortality at age
      SAA_mid[r,y,,] = exp(-0.5 * ZAA[r,y,,]) # Survival at age at midpoint of year

    } # end y loop
  } # end r loop


  ## Recruitment Transformations and Bias Ramp (Methot and Taylor) -------------------------------
  # Mean or virgin recruitment
  if(n_regions > 1) {
    Rec_trans_prop = c(0, Rec_prop) # set up vector for transformation
    Rec_trans_prop = exp(Rec_trans_prop) / sum(exp(Rec_trans_prop)) # do multinomial logit to get recruitment proportions
  } else Rec_trans_prop = 1

  # Global recruitment scalar
  R0 = exp(ln_global_R0) # exponentiate

  # Steepness
  h_trans = 0.2 + (1 - 0.2) * RTMB::plogis(steepness_h) # bound steepness between 0.2 and 1

  # Recruitment SD
  sigmaR2_early = exp(ln_sigmaR[1])^2 # recruitment variability for early period
  sigmaR2_late = exp(ln_sigmaR[2])^2 # recruitment variability for late period

  # Bias ramp set up
  if (do_rec_bias_ramp == 0) {
    bias_ramp = rep(1, n_yrs) # don't do bias ramp, set values to 1
  } else if (do_rec_bias_ramp == 1) {
    # setup bias ramp year ranges
    years = 1:n_yrs # years for indexing
    range1 = which(years >= bias_year[1] & years < bias_year[2])  # ascending limb
    range2 = which(years >= bias_year[2] & years < bias_year[3])  # full bias correction
    range3 = which(years >= bias_year[3] & years < bias_year[4])  # descending limb

    # Apply bias ramp to the different ramp year ranges
    if (length(range1) > 0) bias_ramp[range1] = (years[range1] - bias_year[1]) / (bias_year[2] - bias_year[1]) # ascneding limb
    if (length(range2) > 0) bias_ramp[range2] = 1 # full bias correction
    if (length(range3) > 0) bias_ramp[range3] = 1 - ((years[range3] - bias_year[3]) / (bias_year[4] - bias_year[3])) # descending limb
  } # end if doing bias ramp

  ## Initial Age Structure ---------------------------------------------------
  if(init_age_strc == 0) { # start initial age structure with iterative approach
    for(r in 1:n_regions) {
      for(s in 1:n_sexes) {
        tmp_cumsum_Z = cumsum(natmort[r,1,1:(n_ages-1),s] + init_F * fish_sel[r,1,1:(n_ages-1),s,1])
        Init_NAA[r,,s] = c(R0, R0 * exp(-tmp_cumsum_Z)) * sexratio[s] * Rec_trans_prop[r]
      } # end s loop
    } # end r loop

    # Apply annual cycle and iterate to equilibrium
    for(i in 1:init_iter) {
      for(s in 1:n_sexes) {
        Init_NAA_next_year[,1,s] = R0 * sexratio[s] * Rec_trans_prop # recruitment
        # recruits don't move
        if(do_recruits_move == 0) for(a in 2:n_ages) Init_NAA[,a,s] = t(Init_NAA[,a,s]) %*% Movement[,,1,a,s] # movement
        # recruits move
        if(do_recruits_move == 1) for(a in 1:n_ages) Init_NAA[,a,s] = t(Init_NAA[,a,s]) %*% Movement[,,1,a,s] # movement
        # ageing and mortality
        Init_NAA_next_year[,2:n_ages,s] = Init_NAA[,1:(n_ages-1),s] * exp(-(natmort[,1,1:(n_ages-1),s] + (init_F * fish_sel[,1,1:(n_ages-1),s,1])))
        # accumulate plus group
        Init_NAA_next_year[,n_ages,s] = (Init_NAA_next_year[,n_ages,s] * exp(-(natmort[,1,n_ages,s] + (init_F * fish_sel[,1,n_ages,s,1])))) +
                                        (Init_NAA[,n_ages,s] * exp(-(natmort[,1,n_ages,s] + (init_F * fish_sel[,1,n_ages,s,1]))))
        Init_NAA = Init_NAA_next_year # iterate to next cycle
      } # end s loop
    } # end i loop

    # Apply initial age deviations
    for(r in 1:n_regions) {
      for(s in 1:n_sexes) {
        Init_NAA[r,2:(n_ages-1),s] = Init_NAA[r,2:(n_ages-1),s] * exp(ln_InitDevs[r,]) # add in non-equilibrium age structure
        NAA[r,1,2:n_ages,s] = Init_NAA[r,2:n_ages,s] # add in plus group
      } # end s loop
    } # end r loop
  } # end if

  # Current Assessment Approach -- FLAG, I think it's wrong!
  if(init_age_strc == 1) {
    init_age_idx = 1:(n_ages - 2) # Get initial age indexing
    for(r in 1:n_regions) {
      for(s in 1:n_sexes) {
        NAA[r,1,init_age_idx + 1,s] = R0 * exp(ln_InitDevs[r,init_age_idx] - (init_age_idx * (natmort[r,1, init_age_idx + 1, s] +
                                      (init_F * fish_sel[r,1, init_age_idx + 1, s, 1])))) * sexratio[s] *  Rec_trans_prop[r] # not plus group
        # Plus group calculations
        NAA[r,1,n_ages,s] = R0 * exp( - ((n_ages - 1) * (natmort[r,1, n_ages, s] + (init_F * fish_sel[r,1, n_ages, s, 1]))) ) /
                            (1 - exp(-(natmort[r,1, n_ages, s] + (init_F * fish_sel[r,1, n_ages, s, 1])))) * sexratio[s] * Rec_trans_prop[r]

      } # end s loop
    } # end r loop
  } # end if


  ## Population Projection ---------------------------------------------------
  for(y in 1:n_yrs) {

    ### Annual Recruitment ------------------------------------------------------
    # Get Deterministic Recruitment
    tmp_Det_Rec = Get_Det_Recruitment(recruitment_model = rec_model,
                                      recruitment_dd = rec_dd,
                                      R0 = R0,
                                      Rec_Prop = Rec_trans_prop,
                                      h = h_trans,
                                      n_ages = n_ages,
                                      n_regions = n_regions,
                                      WAA = WAA[,y,,1],
                                      MatAA = MatAA[,y,,1],
                                      natmort = natmort[,y,,1],
                                      SSB_vals = SSB,
                                      y = y,
                                      rec_lag = rec_lag
                                      )
    for(r in 1:n_regions) {
      for(s in 1:n_sexes) {
        if(y < sigmaR_switch) NAA[r,y,1,s] = tmp_Det_Rec[r] * exp(ln_RecDevs[r,y] - (sigmaR2_early/2 * bias_ramp[y])) * sexratio[s] # early period recruitment
        if(y >= sigmaR_switch && y <= n_est_rec_devs) NAA[r,y,1,s] = tmp_Det_Rec[r] * exp(ln_RecDevs[r,y] - (sigmaR2_late/2 * bias_ramp[y])) * sexratio[s] # late period recruitment
        # Dealing with terminal year recruitment
        if(y > n_est_rec_devs) NAA[r,y,1,s] = tmp_Det_Rec[r] * sexratio[s] # mean recruitment in terminal year (not estimate last year rec dev)
      } # end s loop
      Rec[r,y] = sum(NAA[r,y,1,]) # get annual recruitment container here
    } # end r loop

    ### Movement ----------------------------------------------------------------
    if(n_regions > 1) {
      # Recruits don't move
      if(do_recruits_move == 0) {
        # Apply movement after ageing processes - start movement at age 2
        for(a in 2:n_ages) for(s in 1:n_sexes) NAA[,y,a,s] = t(NAA[,y,a,s]) %*% Movement[,,y,a,s]
        for(r in 1:n_regions) NAA[r,y,1,] = Rec[r,y] * sexratio
      } # end if recruits don't move

      # Recruits move here
      if(do_recruits_move == 1) for(a in 1:n_ages) for(s in 1:n_sexes) NAA[,y,a,s] = t(NAA[,y,a,s]) %*% Movement[,,y,a,s]
    } # only compute if spatial

    ### Mortality and Ageing ------------------------------------------------------
    NAA[,y+1,2:n_ages,] = NAA[,y,1:(n_ages-1),] * exp(-ZAA[,y,1:(n_ages-1),]) # Exponential mortality for individuals not in plus group
    NAA[,y+1,n_ages,] = NAA[,y+1,n_ages,] + NAA[,y,n_ages,] * exp(-ZAA[,y,n_ages,]) # Acuumulate plus group

    ### Compute Biomass Quantities ----------------------------------------------
    Total_Biom[,y] = apply(NAA[,y,,,drop = FALSE] * WAA[,y,,,drop = FALSE], 1, sum) # Total biomass
    SSB[,y] = apply(NAA[,y,,1,drop = FALSE] * WAA[,y,,1,drop = FALSE] * MatAA[,y,,1,drop = FALSE], 1, sum) # Spawning Stock Biomass
    if(n_sexes == 1) SSB[,y] = SSB[,y] * 0.5 # If single sex model, multiply SSB calculations by 0.5

  } # end y loop

  ## Fishery Observation Model -----------------------------------------------
  for(r in 1:n_regions) {
    for(y in 1:n_yrs) {
      for(f in 1:n_fish_fleets) {

        fish_q_blk_idx = fish_q_blocks[r,y,f] # get time-block catchability index
        fish_q[r,y,f] = exp(ln_fish_q[r,fish_q_blk_idx,f]) # Input into fishery catchability container

        CAA[r,y,,,f] = FAA[r,y,,,f] / ZAA[r,y,,] * NAA[r,y,,] * (1 - exp(-ZAA[r,y,,])) # Catch at age (Baranov's)

        for(s in 1:n_sexes) {
          if(fit_lengths == 1 && sablefish_ADMB == 1) CAL[r,y,,s,f] = SizeAgeTrans[r,y,,,s] %*% (CAA[r,y,,s,f] / sum(CAA[r,y,,s,f])) # Catch at length (Sablefish bridging specific)
          else if(fit_lengths == 1) CAL[r,y,,s,f] = SizeAgeTrans[r,y,,,s] %*% CAA[r,y,,s,f] # Catch at length
        } # end s loop

        PredCatch[r,y,f] = sum(CAA[r,y,,,f] * WAA[r,y,,]) # get total catch

        # Get fishery index
        if(fish_idx_type[r,f] == 0) PredFishIdx[r,y,f] = fish_q[r,y,f] * sum(NAA[r,y,,] * SAA_mid[r,y,,] * fish_sel[r,y,,,f]) # abundance
        if(fish_idx_type[r,f] == 1) {
          if(fish_q_blk_idx == 1 && sablefish_ADMB == 1) PredFishIdx[r,y,f] = fish_q[r,y,f] * sum(NAA[r,y,,] * SAA_mid[r,y,,] * fish_sel[r,y,,1,f] * WAA[r,y,,]) # first time block (Sablefish bridging specific)
          else PredFishIdx[r,y,f] = fish_q[r,y,f] * sum(NAA[r,y,,] * SAA_mid[r,y,,] * fish_sel[r,y,,,f] * WAA[r,y,,]) # for not first time block
        } # weight

      } # end f loop
    } # end y loop
  } # end r loop

  ## Survey Observation Model ------------------------------------------------
  for(r in 1:n_regions) {
    for(y in 1:n_yrs) {
      for(sf in 1:n_srv_fleets) {

        srv_q_blk_idx = srv_q_blocks[r,y,sf] # get time-block catchability index
        srv_q[r,y,sf] = exp(ln_srv_q[r,srv_q_blk_idx,sf]) # Input into survey catchability container

        if(sablefish_ADMB == 1) SrvIAA[r,y,,,sf] = NAA[r,y,,] * srv_sel[r,y,,,sf] # Survey index at age (sablefish specific)
        else SrvIAA[r,y,,,sf] = NAA[r,y,,] * srv_sel[r,y,,,sf] * SAA_mid[r,y,,] # Survey index at age

        for(s in 1:n_sexes) {
          if(fit_lengths == 1) SrvIAL[r,y,,s,sf] = SizeAgeTrans[r,y,,,s] %*% SrvIAA[r,y,,s,sf] # Survey index at length
        } # end s loop

        # Get predicted survey index
        if(srv_idx_type[r,sf] == 0) PredSrvIdx[r,y,sf] = srv_q[r,y,sf] * sum(NAA[r,y,,] * srv_sel[r,y,,,sf] * SAA_mid[r,y,,]) # abundance
        if(srv_idx_type[r,sf] == 1) PredSrvIdx[r,y,sf] = srv_q[r,y,sf] * sum(NAA[r,y,,] * srv_sel[r,y,,,sf] * SAA_mid[r,y,,] * WAA[r,y,,]) # biomass

      } # end sf loop
    } # end y loop
  } # end r loop


  ## Tagging Observation Model -----------------------------------------------
  if(UseTagging == 1) {

    # Set up tag reporting rates
    for(r in 1:n_regions) {
      TagRep_blk_idx = Tag_Reporting_blocks[r,]  # Get all blocks for this region
      Tag_Reporting[r,] = RTMB::plogis(Tag_Reporting_Pars[r,TagRep_blk_idx])  # inverse logit transform
    } # end r loop

    # Transform tagging parameters here
    Tag_Shed = exp(ln_Tag_Shed) # Transform shedding rates
    Init_Tag_Mort = exp(ln_Init_Tag_Mort) # Transform initial mortality

    # Extract out indexing for tag cohorts
    tr_vec = tag_release_indicator[,1] # tag release region vector
    ty_vec = tag_release_indicator[,2] # tag release year vector

    for(tc in 1:n_tag_cohorts) {
      tr = tr_vec[tc] # extract tag release region
      ty = ty_vec[tc] # extract tag release year

      for(ry in 1:min(max_tag_liberty, n_yrs - ty + 1)) {

        y = ty + ry - 1 # Get index for actual year in the model (instead of tag year)

        # Get total mortality
        tmp_Z = Get_Tagging_Mortality(tag_selex = tag_selex, tag_natmort = tag_natmort,
                                      Fmort = Fmort, natmort = natmort, fish_sel = fish_sel,
                                      Tag_Shed = Tag_Shed, n_regions = n_regions, n_ages = n_ages,
                                      n_sexes = n_sexes, n_fish_fleets = n_fish_fleets, y = y, what = "Z")

        # Get fishing mortality
        tmp_F = Get_Tagging_Mortality(tag_selex = tag_selex, tag_natmort = tag_natmort,
                                      Fmort = Fmort, natmort = natmort, fish_sel = fish_sel,
                                      Tag_Shed = Tag_Shed, n_regions = n_regions, n_ages = n_ages,
                                      n_sexes = n_sexes, n_fish_fleets = n_fish_fleets, y = y, what = "F")

        # Handle tagging dynamics
        if(ry == 1) {
          # Initialize tagging dynamics for first recapture year
          tmp_Z = tmp_Z * t_tagging # discounting mortality if t_tagging != 1
          Tags_Avail[ry,tc,tr,,] = Tagged_Fish[tc,,] * exp(-Init_Tag_Mort) # Input tag releases to the first year
          if(t_tagging == 1) for(a in 1:n_ages) for(s in 1:n_sexes) Tags_Avail[ry,tc,,a,s] = t(Tags_Avail[ry,tc,,a,s]) %*% Movement[,,y,a,s] # Only apply movement if t_tagging == 1 in first recapture year
        } else for(a in 1:n_ages) for(s in 1:n_sexes) Tags_Avail[ry,tc,,a,s] = t(Tags_Avail[ry,tc,,a,s]) %*% Movement[,,y,a,s] # Movement always occurs after first release year

        # Mortality and ageing of tagged fish
        Tags_Avail[ry+1,tc,,2:n_ages,] = Tags_Avail[ry,tc,,1:(n_ages-1),] * exp(-tmp_Z[,1,1:(n_ages-1),]) # not in plus group
        Tags_Avail[ry+1,tc,,n_ages,] = Tags_Avail[ry+1,tc,,n_ages,] + (Tags_Avail[ry,tc,,n_ages,] * exp(-tmp_Z[,1,n_ages,])) # accumulate plus group

        # Get predicted recaptures
        Pred_Tag_Recap[ry,tc,,,] = Tag_Reporting[,y] * (tmp_F[,1,,] / tmp_Z[,1,,]) * Tags_Avail[ry,tc,,,] * (1 - exp(-tmp_Z[,1,,]))

      } # end ry loop
    } # end tc loop
  } # end if for using tagging data


  # Likelihood Equations -------------------------------------------------------------
  ## Fishery Likelihoods -----------------------------------------------------
  ### Fishery Catches ---------------------------------------------------------
  for(y in 1:n_yrs) {
    for(f in 1:n_fish_fleets) {

      # If we have catch data but it's not resolved on the region scale (sum across regions)
      if(UseCatch[1,y,f] == 1 && Catch_Type[y,f] == 0) {
        # ADMB likelihoods
        if(likelihoods == 0) {
          Catch_nLL[1,y,f] = UseCatch[1,y,f] * (log(ObsCatch[1,y,f] + Catch_Constant[f]) -
                                                  log(sum(PredCatch[,y,f]) + Catch_Constant[f]))^2 # SSQ Catch
        } # ADMB likelihoods
        if(likelihoods == 1) {
          Catch_nLL[1,y,f] = UseCatch[1,y,f] -1 * RTMB::dnorm(log(ObsCatch[1,y,f] + Catch_Constant[f]),
                                                              log(sum(PredCatch[,y,f]) + Catch_Constant[f]),
                                                              exp(ln_sigmaC[1,f]), TRUE)
        } # TMB likelihoods
      } # if some fishery catches are aggregated

      for(r in 1:n_regions) {

        # If we have catch data and it's resolved on the region scale
        if(UseCatch[r,y,f] == 1 && Catch_Type[y,f] == 1) {
          # ADMB likelihoods
          if(likelihoods == 0) {
            Catch_nLL[r,y,f] = UseCatch[r,y,f] * (log(ObsCatch[r,y,f] + Catch_Constant[f]) -
                                                    log(PredCatch[r,y,f] + Catch_Constant[f]))^2 # SSQ Catch
          } # ADMB likelihoods
          if(likelihoods == 1) {
            Catch_nLL[r,y,f] = UseCatch[r,y,f] -1 * RTMB::dnorm(log(ObsCatch[r,y,f] + Catch_Constant[f]),
                                                                log(PredCatch[r,y,f] + Catch_Constant[f]),
                                                                exp(ln_sigmaC[r,f]), TRUE)
          } # TMB likelihoods
        } # if no NAs for fishery catches

      } # end r loop
    } # end f loop
  } # end y loop


  ### Fishery Indices ---------------------------------------------------------
  for(y in 1:n_yrs) {
    for(r in 1:n_regions) {
      for(f in 1:n_fish_fleets) {

        # ADMB likelihoods
        if(UseFishIdx[r,y,f] == 1) {
          if(likelihoods == 0) {
            FishIdx_nLL[r,y,f] = (log(ObsFishIdx[r,y,f] + 1e-4) - log(PredFishIdx[r,y,f] + 1e-4))^2 /
              (2 * (ObsFishIdx_SE[r,y,f] / ObsFishIdx[r,y,f])^2) # lognormal Fishery index
          }

          # TMB likelihoods
          if(likelihoods == 1) {
            FishIdx_nLL[r,y,f] = -1 * RTMB::dnorm(log(ObsFishIdx[r,y,f] + 1e-10), log(PredFishIdx[r,y,f] + 1e-10),
                                                  ObsFishIdx_SE[r,y,f], TRUE)
          }
        }

      } # end f loop
    } # end r loop
  } # end y loop

  ### Fishery Compositions ------------------------------------------------
  for(y in 1:n_yrs) {
    for(f in 1:n_fish_fleets) {

      # Fishery Age Compositions
      if(sum(UseFishAgeComps[,y,f]) >= 1) {
        FishAgeComps_nLL[,y,,f] = Get_Comp_Likelihoods(

          # Expected and Observed values
          Exp = CAA[,y,,,f],
          Obs = ObsFishAgeComps[,y,,,f],

          # Input sample size and multinomial weight
          ISS = ISS_FishAgeComps[,y,,f],
          Wt_Mltnml = Wt_FishAgeComps[,y,,f],
          # Composition and Likelihood Type
          Comp_Type = FishAgeComps_Type[y,f],
          Likelihood_Type = FishAgeComps_LikeType[f],

          # overdispersion pars, Number of sexes, regions, age or length comps, and ageing error
          ln_theta = ln_FishAge_theta[,,f],
          ln_theta_agg = ln_FishAge_theta_agg[f],
          LN_corr_pars = FishAge_corr_pars[,,f,],
          LN_corr_pars_agg = FishAge_corr_pars_agg[f],
          n_regions = n_regions, n_sexes = n_sexes, age_or_len = 0,
          AgeingError = AgeingError,
          use = UseFishAgeComps[,y,f],
          n_bins = n_ages,
          comp_agg_type = FishAge_comp_agg_type[f]
        )

      } # if we have fishery age comps

      # Fishery Length Compositions
      if(sum(UseFishLenComps[,y,f]) >= 1 && fit_lengths == 1) {
        FishLenComps_nLL[,y,,f] = Get_Comp_Likelihoods(

          # Expected and Observed values
          Exp = CAL[,y,,,f],
          Obs = ObsFishLenComps[,y,,,f],

          # Input sample size and multinomial weight
          ISS = ISS_FishLenComps[,y,,f],
          Wt_Mltnml = Wt_FishLenComps[,y,,f],

          # Composition and Likelihood Type
          Comp_Type = FishLenComps_Type[y,f],
          Likelihood_Type = FishLenComps_LikeType[f],

          # overdispersion, Number of sexes, regions age or length comps, and ageing error
          ln_theta = ln_FishLen_theta[,,f],
          ln_theta_agg = ln_FishLen_theta_agg[f],
          LN_corr_pars = FishLen_corr_pars[,,f,],
          LN_corr_pars_agg = FishLen_corr_pars_agg[f],
          n_regions = n_regions, n_sexes = n_sexes, age_or_len = 1,
          AgeingError = NA, use = UseFishLenComps[,y,f],
          n_bins = n_lens, comp_agg_type = FishLen_comp_agg_type[f]
        )

      } # if we have fishery length comps

    } # end f loop
  } # end y loop


  ## Survey Likelihoods ------------------------------------------------------
  ### Survey Indices ---------------------------------------------------------
  for(y in 1:n_yrs) {
    for(r in 1:n_regions) {
      for(sf in 1:n_srv_fleets) {

        # ADMB likelihoods
        if(UseSrvIdx[r,y,sf] == 1) {
          if(likelihoods == 0) {
            SrvIdx_nLL[r,y,sf] = (log(ObsSrvIdx[r,y,sf] + 1e-4) - log(PredSrvIdx[r,y,sf] + 1e-4))^2 /
              (2 * (ObsSrvIdx_SE[r,y,sf] / ObsSrvIdx[r,y,sf])^2) # lognormal Survey index
          }

          # TMB likelihoods
          if(likelihoods == 1) {
            SrvIdx_nLL[r,y,sf] = -1 * RTMB::dnorm(log(ObsSrvIdx[r,y,sf] + 1e-10), log(PredSrvIdx[r,y,sf] + 1e-10),
                                                  ObsSrvIdx_SE[r,y,sf], TRUE)
          }
        }

      } # end sf loop
    } # end r loop
  } # end y loop

  ### Survey Compositions ---------------------------------------------------------
  for(y in 1:n_yrs) {
    for(sf in 1:n_srv_fleets) {

      # Survey Age Compositions
      if(sum(UseSrvAgeComps[,y,sf]) >= 1) {
        SrvAgeComps_nLL[,y,,sf] = Get_Comp_Likelihoods(

          # Expected and Observed values
          Exp = SrvIAA[,y,,,sf],
          Obs = ObsSrvAgeComps[,y,,,sf],

          # Input sample size and multinomial weight
          ISS = ISS_SrvAgeComps[,y,,sf],
          Wt_Mltnml = Wt_SrvAgeComps[,y,,sf],

          # Composition and Likelihood Type
          Comp_Type = SrvAgeComps_Type[y,sf],
          Likelihood_Type = SrvAgeComps_LikeType[sf],

          # overdispersion, Number of sexes, regions, age or length comps, and ageing error
          ln_theta = ln_SrvAge_theta[,,sf],
          ln_theta_agg = ln_SrvAge_theta_agg[sf],
          LN_corr_pars = SrvAge_corr_pars[,,sf,],
          LN_corr_pars_agg = SrvAge_corr_pars_agg[sf],
          n_regions = n_regions, n_sexes = n_sexes, age_or_len = 0,
          AgeingError = AgeingError, use = UseSrvAgeComps[,y,sf],
          n_bins = n_ages, comp_agg_type = SrvAge_comp_agg_type[sf]
        )

      } # if we have survey age comps

      # Survey Length Compositions
      if(sum(UseSrvLenComps[,y,sf]) >= 1 && fit_lengths == 1) {
        SrvLenComps_nLL[,y,,sf] = Get_Comp_Likelihoods(

          # Expected and Observed values
          Exp = SrvIAL[,y,,,sf],
          Obs = ObsSrvLenComps[,y,,,sf],

          # Input sample size and multinomial weight
          ISS = ISS_SrvLenComps[,y,,sf],
          Wt_Mltnml = Wt_SrvLenComps[,y,,sf],

          # Composition and Likelihood Type
          Comp_Type = SrvLenComps_Type[y,sf],
          Likelihood_Type = SrvLenComps_LikeType[sf],

          # overdispersion, Number of sexes, regions, age or length comps, and ageing error
          ln_theta = ln_SrvLen_theta[,,sf],
          ln_theta_agg = ln_SrvLen_theta_agg[sf],
          LN_corr_pars = SrvLen_corr_pars[,,sf,],
          LN_corr_pars_agg = SrvLen_corr_pars_agg[sf],
          n_regions = n_regions, n_sexes = n_sexes, age_or_len = 1,
          AgeingError = NA, use = UseSrvLenComps[,y,sf],
          n_bins = n_lens, comp_agg_type = SrvLen_comp_agg_type[sf]
        )

      } # if we have survey length comps

    } # end sf loop
  } # end y loop


  ## Tag Likelihoods ---------------------------------------------------------
  if(UseTagging == 1) {
    for(tc in 1:n_tag_cohorts) {

      # set up tagging cohort indexing
      tr = tr_vec[tc] # extract tag release region
      ty = ty_vec[tc] # extract tag release year

      for(ry in mixing_period:min(max_tag_liberty, n_yrs - ty + 1)) { # loop through recapture years
        for(r in 1:n_regions) {
          for(a in 1:n_move_age_tag_pool) {
            for(s in 1:n_move_sex_tag_pool) {

              move_age_pool_idx = move_age_tag_pool[[a]] # extract movement age pool indices
              move_sex_pool_idx = move_sex_tag_pool[[s]] # extract movement sex pool indices

              # Poisson likelihood
              if(Tag_LikeType == 0) {
                Tag_nLL[ry,tc,r,1,1] = Tag_nLL[ry,tc,r,1,1] + -dpois_noint(sum(Obs_Tag_Recap[ry,tc,r,move_age_pool_idx,move_sex_pool_idx] + 1e-10),
                                                                           sum(Pred_Tag_Recap[ry,tc,r,move_age_pool_idx,move_sex_pool_idx] + 1e-10),
                                                                           give_log = TRUE)
              } # end if poisson likelihood

              # Negative binomial likelihood
              if(Tag_LikeType == 1) {
                log_mu = log(sum(Pred_Tag_Recap[ry,tc,r,move_age_pool_idx,move_sex_pool_idx] + 1e-10)) # log mu
                log_var_minus_mu = 2 * log_mu - ln_tag_theta # log var minus mu
                Tag_nLL[ry,tc,r,1,1] = Tag_nLL[ry,tc,r,1,1] + -dnbinom_robust_noint(x = sum(Obs_Tag_Recap[ry,tc,r,move_age_pool_idx,move_sex_pool_idx] + 1e-10),
                                                                                    log_mu = log_mu, log_var_minus_mu = log_var_minus_mu, give_log = TRUE)
              } # end if for negative binomial likelihood

            } # end s loop
          } # end a loop
        } # end r loop

        # Release Conditioned for Multinomial or Dirichlet-Multinomial
        if(Tag_LikeType %in% c(2, 4)) {

          # Temporary vectors for recaptured individuals
          tmp_pred_c_all = vector()
          tmp_obs_c_all = vector()

          # number of tags released for a given tag cohort
          tmp_n_tags_released = sum(Tagged_Fish[tc,,] + 1e-10)

          # Loop through age and sex pooling and combine vectors into the correct format
          for(a in 1:n_move_age_tag_pool) {
            for(s in 1:n_move_sex_tag_pool) {
              move_age_pool_idx = move_age_tag_pool[[a]] # extract movement age pool indices
              move_sex_pool_idx = move_sex_tag_pool[[s]] # extract movement sex pool indices

              # Pool observed and expected if any pooling
              for (r in 1:n_regions) {
                pred_val = sum(Pred_Tag_Recap[ry, tc, r, move_age_pool_idx, move_sex_pool_idx] + 1e-10) # sum across age and sex groups
                obs_val  = sum(Obs_Tag_Recap[ry, tc, r, move_age_pool_idx, move_sex_pool_idx] + 1e-10) # sum across age and sex groups
                tmp_pred_c_all = c(tmp_pred_c_all, pred_val) # combine predicted recaptures for a given age sex pooled group
                tmp_obs_c_all  = c(tmp_obs_c_all,  obs_val) # combine observed recaptures for a given age sex pooled group
              } # end r loop

            } # end a loop
          } # end s loop

          # Normalize observed and predicted recaptures
          tmp_pred_c_all = tmp_pred_c_all / tmp_n_tags_released
          tmp_obs_c_all = tmp_obs_c_all / tmp_n_tags_released

          # Add in observed and predicted non-recaptures
          tmp_pred = c(tmp_pred_c_all, 1 - sum(tmp_pred_c_all))
          tmp_obs = c(tmp_obs_c_all, 1 - sum(tmp_obs_c_all))
          if(Tag_LikeType == 2) Tag_nLL[ry,tc,1,1,1] = -tmp_n_tags_released * sum((tmp_obs) * log(tmp_pred)) # multinomial
          if(Tag_LikeType == 4) Tag_nLL[ry,tc,1,1,1] =  -1 * ddirmult(obs = tmp_obs, pred = tmp_pred, Ntotal = tmp_n_tags_released, ln_theta = ln_tag_theta, TRUE) # Dirichlet Multinomial
        } # end if release conditioned

        # Recapture Conditioned (Multinomial or Dirichlet-Multinomial)
        if(Tag_LikeType %in% c(3,5)) {
          # Temporary vectors for recaptured individuals
          tmp_pred_all = vector()
          tmp_obs_all = vector()

          # number of recaptures
          tmp_n_tags_recap = sum(Obs_Tag_Recap[ry,tc,,,] + 1e-10)

          # Loop through age and sex pooling and combine vectors into the correct format
          for(a in 1:n_move_age_tag_pool) {
            for(s in 1:n_move_sex_tag_pool) {
              move_age_pool_idx = move_age_tag_pool[[a]] # extract movement age pool indices
              move_sex_pool_idx = move_sex_tag_pool[[s]] # extract movement sex pool indices

              for (r in 1:n_regions) {
                pred_val = sum(Pred_Tag_Recap[ry, tc, r, move_age_pool_idx, move_sex_pool_idx] + 1e-10) # sum across age and sex groups
                obs_val  = sum(Obs_Tag_Recap[ry, tc, r, move_age_pool_idx, move_sex_pool_idx] + 1e-10) # sum across age and sex groups
                tmp_pred_all = c(tmp_pred_all, pred_val) # combine predicted recaptures for a given age sex pooled group
                tmp_obs_all  = c(tmp_obs_all,  obs_val) # combine observed recaptures for a given age sex pooled group
              } # end r loop

            } # end a loop
          } # end s loop

          # Normalize observed and predicted recaptures
          tmp_pred_all = tmp_pred_all / sum(tmp_pred_all)
          tmp_obs_all = tmp_obs_all / tmp_n_tags_recap
          if(Tag_LikeType == 3) Tag_nLL[ry,tc,1,1,1] = -1 * tmp_n_tags_recap * sum(((tmp_obs_all) * log(tmp_pred_all))) # Multinomial
          if(Tag_LikeType == 5) Tag_nLL[ry,tc,1,1,1] =  -1 * ddirmult(obs = tmp_obs_all, pred = tmp_pred_all, Ntotal = tmp_n_tags_recap, ln_theta = ln_tag_theta, TRUE) # Dirichlet Multinomial

        } # end if multinomial recapture conditioned
      } # end ry loop

    } # end tc loop
  } # if we are using tagging data

  ## Priors and Penalties ----------------------------------------------------
  ### Fishing Mortality (Penalty) ---------------------------------------------
  if(Use_F_pen == 1) {
    for(f in 1:n_fish_fleets) {
      for(y in 1:n_yrs) {
        for(r in 1:n_regions) {

          if(UseCatch[r,y,f] == 1) {
            if(likelihoods == 0) {
              if(Catch_Type[y,f] == 0 && est_all_regional_F == 0) Fmort_nLL[1,y,f] = ln_F_devs_AggCatch[y,f]^2 # Use aggregated catch
              else Fmort_nLL[r,y,f] = (ln_F_devs[r,y,f] / exp(ln_sigmaF[r,f]))^2 # SSQ ADMB
            } # ADMB

            if(likelihoods == 1) {
              if(Catch_Type[y,f] == 0 && est_all_regional_F == 0) Fmort_nLL[1,y,f] = -RTMB::dnorm(ln_F_devs_AggCatch[y,f], 0, 1, TRUE) # Use aggregated catch
              else Fmort_nLL[r,y,f] = -RTMB::dnorm(ln_F_devs[r,y,f], 0, exp(ln_sigmaF[r,f]), TRUE)
            } # TMB
          } # end if have catch


        } # end r loop
      } # y loop
    } # f loop
  } #  if using fishing mortality penalty

  ### Selectivity (Penalty) ---------------------------------------------------
  for(r in 1:n_regions) {
    # Fishery Selectivity Deviations
    for(f in 1:n_fish_fleets) {
      if(cont_tv_fish_sel[r,f] > 0) {
        sel_nLL = sel_nLL + - Get_sel_PE_loglik(PE_model = cont_tv_fish_sel[r,f], # process error model
                                                PE_pars = fishsel_pe_pars[,,,f, drop = FALSE], # process error parameters for a given fleet (correlaiton and sigmas)
                                                ln_devs = ln_fishsel_devs[,,,,f, drop = FALSE], # extract out process error deviations for a given fleet
                                                map_sel_devs = map_ln_fishsel_devs[,,,,f, drop = FALSE])
      } # end if

      # Mean Standardizing to help with interpretability
      if(cont_tv_fish_sel[r,f] %in% 3:5) for(s in 1:n_sexes) fish_sel[r,,,s,f] = fish_sel[r,,,s,f] / mean(fish_sel[r,,,s,f])

    } # end f loop

    # Survey Selectivity Deviations
    for(sf in 1:n_srv_fleets) {
      if(cont_tv_srv_sel[r,sf] > 0) {
        sel_nLL = sel_nLL + - Get_sel_PE_loglik(PE_model = cont_tv_srv_sel[r,f], # process error model
                                                PE_pars = srvsel_pe_pars[,,,sf, drop = FALSE], # process error parameters for a given fleet (correlaiton and sigmas)
                                                ln_devs = ln_srvsel_devs[,,,,sf, drop = FALSE], # extract out process error deviations for a given fleet
                                                map_sel_devs = map_ln_srvsel_devs[,,,,sf, drop = FALSE])
      } # end if

      # Mean Standardizing to help with interpretability
      if(cont_tv_srv_sel[r,sf] %in% 3:5) for(s in 1:n_sexes) srv_sel[r,,,s,sf] = srv_sel[r,,,s,sf] / mean(srv_sel[r,,,s,sf])

    } # end sf loop
  } # end r loop

  ### Recruitment (Penalty) ----------------------------------------------------
  if(likelihoods == 0) {
    for(r in 1:n_regions) {
      Init_Rec_nLL[r,] = (ln_InitDevs[r,] / exp(ln_sigmaR[1]))^2 # initial age structure penalty
      if(sigmaR_switch > 1) for(y in 1:(sigmaR_switch-1)) {
        Rec_nLL[r,y] = (ln_RecDevs[r,y]/exp(ln_sigmaR[1]))^2 + bias_ramp[y]*ln_sigmaR[1] # early period
      } # end first y loop
      for(y in (sigmaR_switch:n_est_rec_devs)) {
        Rec_nLL[r,y] = (ln_RecDevs[r,y]/exp(ln_sigmaR[2]))^2 + bias_ramp[y]*ln_sigmaR[2] # late period
      } # end second y loop
    } # end r loop
    Rec_nLL = 0.5 * sum(Rec_nLL)  # multiply by 0.5 and sum
    Init_Rec_nLL = 0.5 * sum(Init_Rec_nLL) # multiply by 0.5 and sum
  }

  if(likelihoods == 1) {
    for(r in 1:n_regions) {
      Init_Rec_nLL[r,] = -RTMB::dnorm(ln_InitDevs[r,], 0, exp(ln_sigmaR[1]), TRUE) # initial age structure penalty
      if(sigmaR_switch > 1) for(y in 1:(sigmaR_switch-1)) {
        Rec_nLL[r,y] = -RTMB::dnorm(ln_RecDevs[r,y], 0, exp(ln_sigmaR[1]), TRUE)
      } # first y loop
      # Note that this penalizes the terminal year rec devs, which is estimated in this case
      for(y in sigmaR_switch:n_est_rec_devs) {
        Rec_nLL[r,y] = -RTMB::dnorm(ln_RecDevs[r,y], 0, exp(ln_sigmaR[2]), TRUE)
      } # end second y loop
    } # end r loop
  }

  ### Fishery Catchability (Prior) -----------------------------------------------
  if(Use_fish_q_prior == 1) {
    unique_fish_q_pars = sort(unique(as.vector(map_fish_q))) # Figure out unique fish q parameters estimated
    for(i in 1:length(unique_fish_q_pars)) {
      par_idx = which(map_fish_q == unique_fish_q_pars[i], arr.ind = TRUE)[1,] # figure out where unique q parameter first occurs
      r = par_idx[1] # get region index
      b = par_idx[2] # get block index
      f = par_idx[3] # get fleet index
      ln_fish_q_val = ln_fish_q[r,b,f] # extract tag reporting rate value
      if(likelihoods == 0) fish_q_nLL = fish_q_nLL + (ln_fish_q_val - log(fish_q_prior[r,b,f,1]))^2 / (2 * (fish_q_prior[r,b,f,2])^2) # ADMB likelihood
      if(likelihoods == 1) fish_q_nLL = fish_q_nLL -RTMB::dnorm(ln_fish_q_val, log(fish_q_prior[r,b,f,1]), fish_q_prior[r,b,f,2], TRUE) # TMB likelihood
    } # end i loop
  }

  ### Survey Catchability (Prior) -----------------------------------------------
  if(Use_srv_q_prior == 1) {
    unique_srv_q_pars = sort(unique(as.vector(map_srv_q))) # Figure out unique srv q parameters estimated
    for(i in 1:length(unique_srv_q_pars)) {
      par_idx = which(map_srv_q == unique_srv_q_pars[i], arr.ind = TRUE)[1,] # figure out where unique q parameter first occurs
      r = par_idx[1] # get region index
      b = par_idx[2] # get block index
      sf = par_idx[3] # get fleet index
      ln_srv_q_val = ln_srv_q[r,b,sf] # extract tag reporting rate value
      if(likelihoods == 0) srv_q_nLL = srv_q_nLL + (ln_srv_q_val - log(srv_q_prior[r,b,sf,1]))^2 / (2 * (srv_q_prior[r,b,sf,2])^2) # ADMB likelihood
      if(likelihoods == 1) srv_q_nLL = srv_q_nLL -RTMB::dnorm(ln_srv_q_val, log(srv_q_prior[r,b,sf,1]), srv_q_prior[r,b,sf,2], TRUE) # TMB likelihood
    } # end i loop
  }

  ### Natural Mortality (Prior) -----------------------------------------------
  if(Use_M_prior == 1) {
    if(likelihoods == 0) M_nLL = (ln_M - log(M_prior[1]))^2 / (2 * (M_prior[2])^2) # ADMB likelihood
    if(likelihoods == 1) M_nLL = -RTMB::dnorm(ln_M, log(M_prior[1]), M_prior[2], TRUE) # TMB likelihood
  } # end if using natural mortality prior

  ### Steepness (Prior) -----------------------------------------------
  if(Use_h_prior == 1) {
    unique_h_pars = sort(unique(as.vector(map_h_Pars))) # Figure out unique steepness parameters estimated
    for(i in 1:length(unique_h_pars)) {
      r = which(map_h_Pars == unique_h_pars[i]) # region index
      tmp_h_beta_pars = get_beta_scaled_pars(low = 0.2, high = 1, mu = h_mu[r], sigma = h_sd[r]) # get alpha and beta parameters
      tmp_h_trans = (h_trans[r] - 0.2) / (1 - 0.2) # transform random variable
      h_nLL = h_nLL - RTMB::dbeta(x = tmp_h_trans, shape1 = tmp_h_beta_pars[1], shape2 = tmp_h_beta_pars[2], log = TRUE) # penalize
    } # end i loop
  } # end if using steepness prior


  ### Movement Rates (Penalty) ------------------------------------------------
  if(cont_vary_movement > 0) {
    Movement_nLL = Movement_nLL + - Get_move_PE_loglik(PE_model = cont_vary_movement,
                                                       PE_pars = move_pe_pars,
                                                       logit_devs = logit_move_devs,
                                                       map_move_devs = map_logit_move_devs,
                                                       do_recruits_move = do_recruits_move
                                                       )
  }

  ### Movement Rates (Prior) ------------------------------------------------
  # NOTE: If continuous varying movement is estimated, there should only be one set of movement parameters
  # estimated (i.e., the base, mean movement parameters), such that the prior is applied onto the base parameters
  if(Use_Movement_Prior == 1) {
    unique_movement_pars = sort(unique(as.vector(map_Movement_Pars))) # Figure out unique movement parameters estimated
    for(i in 1:length(unique_movement_pars)) {
      par_idx = which(map_Movement_Pars == unique_movement_pars[i], arr.ind = TRUE)[1,] # figure out where unique movement parameter first occurs
      r_from = par_idx[1] # from region
      y = par_idx[3] # year index
      a = par_idx[4] # age index
      s = par_idx[5] # sex index
      Movement_nLL = Movement_nLL - ddirichlet(x = Movement[r_from,,y,a,s], alpha = Movement_prior[r_from,,y,a,s], log = TRUE) # dirichlet prior
    } # end i loop
  } # end if using movement prior


  ### Tag Reporting Rate (Prior) --------------------------------------------
  if(Use_TagRep_Prior == 1) {
    unique_tagrep_pars = sort(unique(as.vector(map_Tag_Reporting_Pars))) # Figure out unique tag reporting parameters estimated
    for(i in 1:length(unique_tagrep_pars)) {
      par_idx = which(map_Tag_Reporting_Pars == unique_tagrep_pars[i], arr.ind = TRUE)[1,] # figure out where unique tagrep parameter first occurs
      r = par_idx[1] # get region index
      b = par_idx[2] # get block index
      TagRep_val = RTMB::plogis(Tag_Reporting_Pars[r,b]) # extract tag reporting rate value
      if(TagRep_PriorType == 0) {
        TagRep_nLL = TagRep_nLL - dbeta_symmetric(p_val = TagRep_val, p_ub = 1, p_lb = 0, p_prsd = TagRep_sd, log = TRUE) # penalize
      } # end if symmetric beta
      if(TagRep_PriorType == 1) {
        a = TagRep_mu / (TagRep_sd * TagRep_sd) # alpha parameter
        b = (1 - TagRep_mu) / (TagRep_sd * TagRep_sd) # beta parameter
        TagRep_nLL = TagRep_nLL -RTMB::dbeta(x = TagRep_val, shape1 = a, shape2 = b, log = TRUE) # penalize
      } # end if for full beta
    } # end i loop
  } # if use tag reporting prior

  # Apply likelihood weights here and compute joint negative log likelihood
  jnLL = (Wt_Catch * sum(Catch_nLL)) + # Catch likelihoods
    (Wt_FishIdx * sum(FishIdx_nLL)) + # Fishery Index likelihood
    (Wt_SrvIdx * sum(SrvIdx_nLL)) + # Survey Index likelihood
    sum(FishAgeComps_nLL) + # Fishery Age likelihood
    sum(FishLenComps_nLL) + # Fishery Length likelihood
    sum(SrvAgeComps_nLL) + # Survey Age likelihood
    sum(SrvLenComps_nLL) + # Survey Length likelihood
    (Wt_Tagging * sum(Tag_nLL)) + # Tagging likelihood
    (Wt_F * sum(Fmort_nLL)) + # Fishery Mortality Penalty
    (Wt_Rec * sum(Rec_nLL)) + # Recruitment Penalty
    (Wt_Rec * sum(Init_Rec_nLL)) + #  Initial Age Penalty
    sel_nLL + #  selectivity penalty
    M_nLL + # Natural Mortality Prior
    h_nLL + # Steepness Prior
    Movement_nLL + # movement Prior
    TagRep_nLL + # tag reporting rate Prior
    fish_q_nLL + # fishery q prior
    srv_q_nLL  # survey q prior

  # Report Section ----------------------------------------------------------
  # Biological Processes
  RTMB::REPORT(R0)
  RTMB::REPORT(h_trans)
  RTMB::REPORT(NAA)
  RTMB::REPORT(ZAA)
  RTMB::REPORT(natmort)
  RTMB::REPORT(bias_ramp)
  RTMB::REPORT(Movement)

  # Fishery Processes
  RTMB::REPORT(init_F)
  RTMB::REPORT(Fmort)
  RTMB::REPORT(FAA)
  RTMB::REPORT(CAA)
  RTMB::REPORT(CAL)
  RTMB::REPORT(PredCatch)
  RTMB::REPORT(PredFishIdx)
  RTMB::REPORT(fish_sel)
  RTMB::REPORT(fish_q)

  # Survey Processes
  RTMB::REPORT(PredSrvIdx)
  RTMB::REPORT(srv_sel)
  RTMB::REPORT(srv_q)
  RTMB::REPORT(SrvIAA)
  RTMB::REPORT(SrvIAL)

  # Tagging Processes
  if(UseTagging == 1) {
    RTMB::REPORT(Pred_Tag_Recap)
    RTMB::REPORT(Tags_Avail)
    RTMB::REPORT(Tag_Reporting)
  }

  # Likelihoods
  RTMB::REPORT(Catch_nLL)
  RTMB::REPORT(FishIdx_nLL)
  RTMB::REPORT(SrvIdx_nLL)
  RTMB::REPORT(FishAgeComps_nLL)
  RTMB::REPORT(SrvAgeComps_nLL)
  RTMB::REPORT(FishLenComps_nLL)
  RTMB::REPORT(SrvLenComps_nLL)
  RTMB::REPORT(M_nLL)
  RTMB::REPORT(Fmort_nLL)
  RTMB::REPORT(Rec_nLL)
  RTMB::REPORT(Init_Rec_nLL)
  RTMB::REPORT(Rec_nLL)
  RTMB::REPORT(Tag_nLL)
  RTMB::REPORT(h_nLL)
  RTMB::REPORT(fish_q_nLL)
  RTMB::REPORT(sel_nLL)
  RTMB::REPORT(srv_q_nLL)
  RTMB::REPORT(Movement_nLL)
  RTMB::REPORT(TagRep_nLL)
  RTMB::REPORT(jnLL)

  # Report for derived quantities
  RTMB::REPORT(Total_Biom)
  RTMB::REPORT(SSB)
  RTMB::REPORT(Rec)

  # Report these in log space because can't be < 0
  RTMB::ADREPORT(log(Total_Biom))
  RTMB::ADREPORT(log(SSB))
  RTMB::ADREPORT(log(Rec))
  RTMB::ADREPORT(Total_Biom)
  RTMB::ADREPORT(SSB)
  RTMB::ADREPORT(Rec)

  return(jnLL)
} # end function
