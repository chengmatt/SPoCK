#' Simulates a spatial, sex, and age-structured population and takes a variety of inputs, with dependecies on various functions (see example)
#' @param output_path path to output simulation objects
#'
#' @returns a list object with a bunch of simulated values and outputs
#' @export Simulate_Pop
#' @import compResidual
#' @import stats
#' @examples
#' \dontrun{
#' library(here)
#' # Setup Simulation Parameters --------------------------------------------
#' # Set up model dimensions
#' Setup_Sim_Dim(n_sims = 100,
#'               n_yrs = 30,
#'               n_regions = 2,
#'               n_ages = 10,
#'               n_sexes = 1,
#'               n_fish_fleets = 1,
#'               n_srv_fleets = 1)
#'
#' Setup_Sim_Containers(n_sims = n_sims,
#'                      n_yrs = n_yrs,
#'                      n_regions = n_regions,
#'                      n_ages = n_ages,
#'                      n_sexes = n_sexes,
#'                      n_fish_fleets = n_fish_fleets,
#'                      n_srv_fleets = n_srv_fleets)
#'
#' # Setup fishing mortality
#' Setup_Sim_FishMort(n_sims = n_sims,
#'                    n_yrs = n_yrs,
#'                    n_regions = n_regions,
#'                    n_fish_fleets = n_fish_fleets,
#'                    sigmaC = 1e-3,
#'                    init_F = matrix(0, nrow = n_regions, ncol = n_fish_fleets),
#'                    Fmort_pattern = matrix(c('constant', "two-way"), nrow = n_regions, ncol = n_fish_fleets),
#'                    Fmort_start = matrix(c(0.5, 0.5), nrow = n_regions, ncol = n_fish_fleets),
#'                    Fmort_fct = matrix(c(1, 0.5), nrow = n_regions, ncol = n_fish_fleets),
#'                    proc_error = TRUE,
#'                    proc_error_sd = 0.1)
#'
#' # Setup fishery selectivity
#' Setup_Sim_FishSel(n_sims = n_sims,
#'                   n_yrs = n_yrs,
#'                   n_regions = n_regions,
#'                   n_ages = n_ages,
#'                   n_sexes = n_sexes,
#'                   n_fish_fleets = n_fish_fleets,
#'                   sel_model = matrix(c('logistic', "logistic"), nrow = n_regions, ncol = n_fish_fleets),
#'                   # a50, k for logistic shared across regions
#'                   fixed_fish_sel_pars = array(c(5,5,1,1), dim = c(n_regions, n_sexes, n_fish_fleets, 2))
#' )
#'
#' # Setup survey catchability and selectivity
#' Setup_Sim_Survey(n_sims = n_sims,
#'                  n_yrs = n_yrs,
#'                  n_regions = n_regions,
#'                  n_ages = n_ages,
#'                  n_sexes = n_sexes,
#'                  n_srv_fleets = n_srv_fleets,
#'                  sigmaSrvIdx = array(0.2, dim = c(n_regions, n_srv_fleets)), # survey observation error
#'                  base_srv_q = array(1, dim = c(n_regions, n_srv_fleets)), # base survey catchability value
#'                  srv_q_pattern = matrix(c('constant', "constant"),
#'                                         nrow = n_regions, ncol = n_srv_fleets), # catchability pattern
#'                  sel_model = matrix(c('logistic', "logistic"),
#'                                     nrow = n_regions, ncol = n_srv_fleets), # selectivity model
#'                  # a50, k, for logistic shared across regions
#'                  fixed_srv_sel_pars = array(c(5,5,1,1), dim = c(n_regions, n_sexes, n_srv_fleets, 2))
#' )
#'
#' # Setup recruitment stuff
#' Setup_Sim_Rec(n_yrs = n_yrs,
#'               n_regions = n_regions,
#'               n_sexes = n_sexes,
#'               n_sims = n_sims,
#'               do_recruits_move = 0, # == 0, recruits don't move , == 1 recruits move
#'               base_rec_sexratio = 1, # single sex
#'               rec_sexratio_vary = "constant",
#'               base_r0 = c(50, 10),
#'               r0_vary = "constant",
#'               base_h = 1,
#'               init_sigmaR = 0.5,
#'               sigmaR = 0.5,
#'               recruitment_opt = 0, # == 0, mean recruitment, == 1 BH
#'               recdev_opt = 0) # == 0, global rec devs, == 1, local rec devs
#'
#' # Setup biologicals
#' Setup_Sim_Biologicals(n_sims = n_sims,
#'                       n_yrs = n_yrs,
#'                       n_regions = n_regions,
#'                       n_ages = n_ages,
#'                       n_sexes = n_sexes,
#'                       base_M_value = array(0.15, dim = c(n_regions, n_ages, n_sims)),
#'                       M_pattern = "constant",
#'                       base_WAA_values = array(rep(5 * (1 - exp(-0.1 * 1:n_ages)),
#'                                                   each = n_regions * n_sexes),
#'                                               dim = c(n_regions, n_ages, n_sexes)),
#'                       WAA_pattern = "constant",
#'                       base_Maturity_AA_values = array(rep(1 / (1 + exp(-0.3 * 1:n_ages)),
#'                                                           each = n_regions * n_sexes),
#'                                                       dim = c(n_regions, n_ages, n_sexes)),
#'                       Maturity_AA_pattern = "constant"
#' )
#'
#' # Setup tagging stuff
#' Setup_Sim_Tagging(n_sims = n_sims,
#'                   n_yrs = n_yrs,
#'                   n_regions = n_regions,
#'                   n_ages = n_ages,
#'                   n_sexes = n_sexes,
#'                   n_tags = 1e4,
#'                   max_liberty = 30,
#'                   tag_years = seq(1, n_yrs, 3),
#'                   t_tagging = 0.5,
#'                   base_Tag_Reporting = rep(0.2, n_regions),
#'                   Tag_Reporting_pattern = "constant",
#'                   Tag_Ind_Mort = 0,
#'                   Tag_Shed = 0
#' )
#'
#' # Setup observation processes
#' Setup_Sim_Observation_Proc(n_fish_fleets = n_fish_fleets,
#'                            n_srv_fleets = n_srv_fleets,
#'                            Comp_Structure = "JntSex_JntRegion",
#'                            Comp_Srv_Like = "Multinomial",
#'                            Comp_Fish_Like = "Multinomial",
#'                            Srv_Like_Pars = NA,
#'                            Fish_Like_Pars = NA,
#'                            Tag_Like = "Poisson",
#'                            Tag_Like_Pars = NA
#' )
#'
#' # Setup movement matrix
#' movement_matrix <- array(0, dim = c(n_regions, n_regions, n_yrs, n_ages, n_sexes, n_sims)) # From, To
#' for(y in 1:n_yrs) {
#'   for(a in 1:n_ages) {
#'     for(s in 1:n_sexes) {
#'       for(sim in 1:n_sims) {
#'         movement_matrix[,,y,,s,sim] <- rep(0.5, 4)
#'       }
#'     } # end s loop
#'   } # end a loop
#' } # end y loop
#'
#' # Run Simulation ----------------------------------------------------------
#' Simulate_Pop(output_path = here("sim_out.RDS"))
#' }
Simulate_Pop <- function(output_path) {

  for(sim in 1:n_sims) {
    for(y in 1:n_yrs) {

      # Initialize Age Structure ------------------------------------------------
      if(y == 1) {

        # Set up initial equilibrium age structure, with cumulative sum of selectivity incorporated
        for(r in 1:n_regions) {
          for(s in 1:n_sexes) {
            tmp_cumsum_Z = cumsum(M[1,r,1:(n_ages-1),s,sim] + init_F[1,r,1,sim] * fish_sel[1,r,1:(n_ages-1),s,1,sim])
            Init_NAA[r,,s,sim] = c(r0[1,r,sim], r0[1,r,sim] * exp(-tmp_cumsum_Z)) * rec_sexratio[1,r,s,sim]
          } # end s loop
        } # end r loop

        # Apply annual cycle and iterate to equilibrium
        for(i in 1:init_iter) {
          for(s in 1:n_sexes) {
            Init_NAA_next_year[,1,s,sim] = r0[1,,sim] * rec_sexratio[1,,s,sim] # recruitment

            # recruits don't move
            if(do_recruits_move == 0) for(a in 2:n_ages) Init_NAA[,a,s,sim] = t(Init_NAA[,a,s,sim]) %*% movement_matrix[,,1,a,s,sim] # movement
            # recruits move
            if(do_recruits_move == 1) for(a in 1:n_ages) Init_NAA[,a,s,sim] = t(Init_NAA[,a,s,sim]) %*% movement_matrix[,,1,a,s,sim] # movement

            # ageing and mortality
            Init_NAA_next_year[,2:n_ages,s,sim] = Init_NAA[,1:(n_ages-1),s,sim] * exp(-(M[1,,1:(n_ages-1),s,sim] + (init_F[1,,1,sim] * fish_sel[1,,1:(n_ages-1),s,1,sim])))
            # accumulate plus group
            Init_NAA_next_year[,n_ages,s,sim] = (Init_NAA_next_year[,n_ages,s,sim] * exp(-(M[1,,n_ages,s,sim] + (init_F[1,,1,sim] * fish_sel[1,,n_ages,s,1,sim])))) +
              (Init_NAA[,n_ages,s,sim] * exp(-(M[1,,n_ages,s,sim] + (init_F[1,,1,sim] * fish_sel[1,,n_ages,s,1,sim]))))
            Init_NAA = Init_NAA_next_year # iterate to next cycle
          } # end s loop
        } # end i loop

        # Apply initial age structure deviations here (FLAG: Revise to incorporate more options)
        tmp_ln_init_devs <- NULL # Initialize container vector to allow for global recruitment
        for(r in 1:n_regions) {
          if(recdev_opt == 0 && is.null(tmp_ln_init_devs)) tmp_ln_init_devs <- stats::rnorm(n_ages-2, 0, init_sigmaR[r]) # simulate initial deviations (global density dependence)
          if(recdev_opt == 1) tmp_ln_init_devs <- stats::rnorm(n_ages-2, 0, init_sigmaR[r]) # simulate initial deviations (local density dependence)
          Init_NAA[r,2:(n_ages-1),s,sim] <- Init_NAA[r,2:(n_ages-1),s,sim] * rep(exp(tmp_ln_init_devs - init_sigmaR[r]^2/2), n_sexes) # apply deviations
          NAA[1,r,2:n_ages,,sim] <- Init_NAA[r,2:n_ages,s,sim] # Plug in initial age structure into 1st year (w/o recruitment)
        } # end r loop
      } # end initializing age structure

      # Apply Movement ----------------------------------------------------------
      for(a in 1:n_ages) for(s in 1:n_sexes) NAA[y,,a,s,sim] <- NAA[y,,a,s,sim] %*% movement_matrix[,,y,a,s,sim]

      # Run Annual Cycle --------------------------------------------------------
      tmp_ln_rec_devs <- NULL # Initialize container vector to allow for global recruitment (remains NULL within a given year)
      for(r in 1:n_regions) {

        ### Recruitment -------------------------------------------------------------
        if(recdev_opt == 0 && is.null(tmp_ln_rec_devs)) tmp_ln_rec_devs <- stats::rnorm(1, 0, sigmaR[y,r]) # Get recruitment deviates (global density dependence)
        if(recdev_opt == 1) tmp_ln_rec_devs <- ln_rec_devs[y,r,sim] <- stats::rnorm(1, 0, sigmaR[y,r]) # Get recruitment deviates (local density dependence)
        ln_rec_devs[y,r,sim] = tmp_ln_rec_devs # input vector of temporary rec devs

        # get deterministic recruitment
        tmp_Det_Rec <- Get_Det_Recruitment(recruitment_model = recruitment_opt, R0 = r0[y,r,sim], h = h[y,r,sim], n_ages = n_ages,
                                           WAA = WAA[y,r,,1,sim], MatAA = Maturity_AA[y,r,,1,sim], natmort = M[y,r,,1,sim],
                                           SSB_vals = SSB[y,r,sim], y = y, rec_lag = rec_lag)

        for(s in 1:n_sexes) {
          # get recruitment
          NAA[y,r,1,s,sim] <- r0[y,r,sim] * exp(ln_rec_devs[y,r,sim] - sigmaR[y,r]^2/2) * rec_sexratio[y,r,s,sim]

          ### Movement ----------------------------------------------------------------
          # Recruits don't move
          if(do_recruits_move == 0) for(a in 2:n_ages) for(s in 1:n_sexes) NAA[y,,a,s,sim] = t(NAA[y,,a,s,sim]) %*% movement_matrix[,,y,a,s,sim]
          # Recruits move here
          if(do_recruits_move == 1) for(a in 1:n_ages) for(s in 1:n_sexes) NAA[y,,a,s,sim] = t(NAA[y,,a,s,sim]) %*% movement_matrix[,,y,a,s,sim]

          ### Mortality and Ageing ----------------------------------------------------
          for(a in 1:n_ages) {
            # Get total mortality
            Z[y,r,a,s,sim] <- M[y,r,a,s,sim] + sum(Fmort[y,r,,sim] * fish_sel[y,r,a,s,,sim]) # Z = M + Fmort
            if(a < n_ages) {
              # Exponential mortality for individuals not in plus group (recruits experience mortality )
              NAA[y+1,r,a+1,s,sim] <- NAA[y,r,a,s,sim] * exp(-Z[y,r,a,s,sim])
            } else {
              # Accumulate individuals recently "recruited" into plus group and individuals from previous year
              NAA[y+1,r,n_ages,s,sim] <- NAA[y+1,r,n_ages,s,sim] + NAA[y,r,n_ages,s,sim] * exp(-Z[y,r,a,s,sim])
            } # end else (calculations for plus group)
          } # end a loop
        } # end s loop

        ### Compute Biomass Quantities ----------------------------------------------
        Total_Biom[y,r,sim] <- sum(as.vector(NAA[y,r,,,sim]) * as.vector(WAA[y,r,,,sim])) # Total Biomass
        SSB[y,r,sim] <- sum(as.vector(NAA[y,r,,1,sim]) * as.vector(WAA[y,r,,1,sim]) * Maturity_AA[y,r,,1,sim]) # Spawning Stock Biomass
        if(n_sexes == 1) SSB[y,r,sim] <- SSB[y,r,sim] * 0.5 # If single sex model, multiply SSB calculations by 0.5

      } # end r loop

    } # end y loop

    # Catch Accounting --------------------------------------------------------
    for(y in 1:n_yrs) {
      for(f in 1:n_fish_fleets) {
        for(r in 1:n_regions) {
          for(s in 1:n_sexes) {
            # Baranov's catch equation
            CAA[y,r,,s,f,sim] <- (Fmort[y,r,f,sim] * fish_sel[y,r,,s,f,sim]) / Z[y,r,,s,sim] *  NAA[y,r,,s,sim] * (1 - exp(-Z[y,r,,s,sim]))

            # Structuring composition data to be split by region and sex
            if(comp_strc == 0) {
              tmp_FishAgeComps_Prob <- CAA[y,r,,s,f,sim] # Get probabilities for a given region and sex
              if(comp_fish_like[f] == 0) Obs_FishAgeComps[y,r,,s,f,sim] <- array(stats::rmultinom(1, 400, tmp_FishAgeComps_Prob / sum(tmp_FishAgeComps_Prob)), dim = c(dim(CAA[y,r,,s,f,sim]))) # simulate multinomial probabilities
              if(comp_fish_like[f] == 1) Obs_FishAgeComps[y,r,,s,f,sim] <- array(compResidual::rdirM(1, 400, 5 * 400 * tmp_FishAgeComps_Prob / sum(tmp_FishAgeComps_Prob)), dim = dim(tmp_FishAgeComps_Prob)) # Simulate dirichlet multinomial samples
            } # end if for "Split" approach for composition data (split by region and sex

          } # end s loop
          # Generate Catch Data
          True_Catch[y,r,f,sim] <- sum(CAA[y,r,,,f,sim] * WAA[y,r,,,sim]) # True Catch
          Obs_Catch[y,r,f,sim] <- True_Catch[y,r,f,sim] * exp(stats::rnorm(1, 0, sigmaC)) # Observed Catch w/ lognormal deviations

          # Structuring composition data to be split by reigon
          if(comp_strc == 1) {
            tmp_FishAgeComps_Prob <- CAA[y,r,,,f,sim] # Get probabilities for a given region and sex
            if(comp_fish_like[f] == 0) Obs_FishAgeComps[y,r,,,f,sim] <- array(stats::rmultinom(1, 400, tmp_FishAgeComps_Prob / sum(tmp_FishAgeComps_Prob)), dim = c(dim(CAA[y,r,,,f,sim, drop = FALSE]))) # simulate multinomial probabilities
            if(comp_fish_like[f] == 1) Obs_FishAgeComps[y,r,,,f,sim] <- array(compResidual::rdirM(1, 400, 5 * 400 * tmp_FishAgeComps_Prob / sum(tmp_FishAgeComps_Prob)), dim = dim(tmp_FishAgeComps_Prob)) # Simulate dirichlet multinomial samples
          } # end if for 'Split' approach for composition data split by region by not by sex

        } # end r loop

        # Structuring composition data to be joint across regions, ages, and sexes
        if(comp_strc == 2) {
          # Store temporary probabilities ordered by ages, sexes, regions (i.e., age 1-30, sex 1, region 1, age 1-30, sex 2, region 1,
          # age 1-30, sex 1, region 2, age 1-30, sex 2, region 2 ... )
          tmp_FishAgeComps_Prob <- aperm(CAA[y, , , , f, sim, drop = FALSE], c(3,4,2,1,5,6)) # ordered by ages, sexes, regions, year = y, fishery fleet = f, and sim = sim
          if(comp_fish_like[f] == 0) tmp_sim_FishAgeComps <- array(stats::rmultinom(1, 500, tmp_FishAgeComps_Prob / sum(tmp_FishAgeComps_Prob)), dim = dim(tmp_FishAgeComps_Prob)) # Simulate Multinomial samples
          if(comp_fish_like[f] == 1) tmp_sim_FishAgeComps <- array(compResidual::rdirM(1, 500, 1 * 500 * tmp_FishAgeComps_Prob / sum(tmp_FishAgeComps_Prob)), dim = dim(tmp_FishAgeComps_Prob)) # Simulate dirichlet multinomial samples
          if(comp_fish_like[f] == 2) tmp_sim_FishAgeComps <- array(rlogistnormal(exp = tmp_FishAgeComps_Prob, pars = 3, comp_like = comp_fish_like[f]), dim = dim(tmp_FishAgeComps_Prob)) # Simulate logistic normal samples
          if(comp_fish_like[f] == 5) tmp_sim_FishAgeComps <- array(rlogistnormal(exp = tmp_FishAgeComps_Prob, pars = c(0.3, 0.6, 0.6, 0.6), comp_like = comp_fish_like[f]), dim = dim(tmp_FishAgeComps_Prob)) # Simulate logistic normal samples
          # Inputing simulated data into dataframe while reshaping to correct dimension (revert to year, region, ages, sexes, fleet, sim)
          Obs_FishAgeComps[y,,,,f,sim] <- aperm(tmp_sim_FishAgeComps, c(4,3,1,2,5,6))
        } # end if for "Joint" approach for composition data across regions, ages, and sexes

      } # end f loop
    } # end y loop

    # Survey Accounting -------------------------------------------------------
    for(y in 1:n_yrs) {
      for(sf in 1:n_srv_fleets) {
        for(r in 1:n_regions) {
          for(s in 1:n_sexes) {
            # Survey Ages Indexed
            Srv_IAA[y,r,,s,sf,sim] <- NAA[y,r,,s,sim] * srv_sel[y,r,,s,sf,sim] * exp(-0.5 * Z[y,r,,s,sim])

            # Structuring composition data to be split by region and sex
            if(comp_strc == 0) {
              tmp_SrvAgeComps_Prob <- Srv_IAA[y,r,,s,sf,sim] # Get probabilities for a given region and sex
              if(comp_srv_like[sf] == 0) Obs_SrvAgeComps[y,r,,s,sf,sim] <- array(stats::rmultinom(1, 400, tmp_SrvAgeComps_Prob / sum(tmp_SrvAgeComps_Prob)), dim = c(dim(Srv_IAA[y,r,,s,sf,sim]))) # simulate multinomial probabilities
            } # end if for "Split" approach for composition data (split by region and sex)

          } # end s loop
          True_SrvIdx[y,r,sf,sim] <- srv_q[y,r,sf,sim] * sum(Srv_IAA[y,r,,,sf,sim]) # True Survey Index
          Obs_SrvIdx[y,r,sf,sim] <- srv_q[y,r,sf,sim] * True_SrvIdx[y,r,sf,sim] * exp(stats::rnorm(1, 0, sigmaSrvIdx[r,sf])) # Observed survey index w/ lognormal deviations

          # Structuring composition data to be split by reigon
          if(comp_strc == 1) {
            tmp_SrvAgeComps_Prob <- Srv_IAA[y,r,,,sf,sim] # Get probabilities for a given region and sex
            if(comp_srv_like[sf] == 0) Obs_SrvAgeComps[y,r,,,sf,sim] <- array(stats::rmultinom(1, 400, tmp_SrvAgeComps_Prob / sum(tmp_SrvAgeComps_Prob)), dim = c(dim(Srv_IAA[y,r,,,sf,sim,drop=FALSE]))) # simulate multinomial probabilities
          } # end if for 'Split' approach for composition data split by region by not by sex
        } # end r loop

        # Structuring composition data to be joint across regions, ages, and sexes
        if(comp_strc == 2) {
          # Store temporary Probabilities ordered by ages, sexes, regions
          # (i.e., age 1-30, sex 1, region 1, age 1-30, sex 2, region 1, age 1-30, sex 1, region 2, age 1-30, sex 2, region 2 ... )
          tmp_SrvAgeComps_Prob <- aperm(Srv_IAA[y, , , , sf, sim, drop = FALSE], c(3,4,2,1,5,6)) # ordered by ages, sexes, regions, year = y, survey fleet = sf, and sim = sim
          if(comp_srv_like[sf] == 0) tmp_sim_SrvAgeComps <- array(stats::rmultinom(1, 500, tmp_SrvAgeComps_Prob / sum(tmp_SrvAgeComps_Prob)), dim = dim(tmp_SrvAgeComps_Prob)) # Simulate Multinomial samples
          # if(comp_srv_like == 1) tmp_sim_SrvAgeComps <- array(compResidual::rdirM(1, 400, 400 * 1 * tmp_SrvAgeComps_Prob / sum(tmp_SrvAgeComps_Prob)), dim = dim(tmp_SrvAgeComps_Prob)) # Simulate Multinomial samples
          # Inputing simulated data into dataframe while reshaping to correct dimension (revert to year, region, ages, sexes, fleet, sim)
          Obs_SrvAgeComps[y,,,,sf,sim] <- aperm(tmp_sim_SrvAgeComps, c(4,3,1,2,5,6))
        } # end if for "Joint" approach for composition data across regions, ages, and sexes

      } # end sf loop
    } # end y loop

    # Tag Releases and Recoveries ------------------------------------------------------------
    for(tag_rel in 1:n_tag_rel_events) {

      # Tag Releases
      tag_rel_yr <- tag_rel_indicator$tag_yrs[tag_rel] # get tag release year
      tag_rel_region <- tag_rel_indicator$regions[tag_rel] # get tag release region
      n_tags_rel <- round(Obs_SrvIdx[tag_rel_yr,,1,sim] / sum(Obs_SrvIdx[tag_rel_yr,,1,sim]) * n_tags) # distribute tags relative to regional survey abundance
      tmp_SrvAgeComps_Prob <- as.vector(Srv_IAA[tag_rel_yr, tag_rel_region, , , 1, sim]) # Use survey proportions to distribute tags
      tagged_fish <- round((tmp_SrvAgeComps_Prob / sum(tmp_SrvAgeComps_Prob)) * n_tags_rel[tag_rel_region]) # Distribute tags across ages, sexes, for a given release event
      Tag_Fish[tag_rel,,,sim] <- array(tagged_fish, dim = c(n_ages, n_sexes)) # Reshape format by ages, and sexes and input into other array

      # Tag Recaptures (do we need to pool individuals into a plus max liberty group or can we just ignore?)
      for(recap_yr in 1:min(max_liberty, n_yrs - tag_rel_yr + 1)) { # recapture year for a given cohort - adding a cut off for time at liberty so we are not tracking fish for 400+ years
        actual_yr <- tag_rel_yr + recap_yr - 1 # Define actual year for indexing purposes

        # Input tagged fish into available tags for recapture and adjust initial number of tagged fish for tag induced mortality (exponential mortality process)
        if(recap_yr == 1) Tag_Avail[1,tag_rel,tag_rel_region,,,sim] <- Tag_Fish[tag_rel,,,sim] * exp(-Tag_Ind_Mort)

        # Get mortality estimates and account for tag shedding
        tmp_F <- Fmort[actual_yr,,1,sim] # temporary fishing mortality variable (uniform sel)
        if(recap_yr == 1 && t_tagging != 1) tmp_Z <- (M[actual_yr,,,,sim, drop = FALSE] + tmp_F + Tag_Shed) * t_tagging  # temporary total mortality variable
        else tmp_Z <- (M[actual_yr,,,,sim, drop = FALSE] + tmp_F + Tag_Shed)

        # Move tagged fish around (movement only occurs after first recapture year if tagging happens midyear, since movement happens at start of yr)
        if(t_tagging != 1 && recap_yr == 1) {} else for(a in 1:n_ages) for(s in 1:n_sexes) Tag_Avail[recap_yr,tag_rel,,a,s,sim] <- Tag_Avail[recap_yr,tag_rel,,a,s,sim] %*% movement_matrix[,,actual_yr,a,s,sim] # recruits move

        # Apply mortality and ageing to tagged fish
        for(a in 1:n_ages) {
          for(s in 1:n_sexes) {
            if(a < n_ages) { # If not in plus group
              # same dynamics as all other fish, but with uniform selex
              Tag_Avail[recap_yr+1,tag_rel,,a+1,s,sim] <- Tag_Avail[recap_yr,tag_rel,,a,s,sim]  * exp(-tmp_Z[1,,a,s,1])
            } else{ # Accumulate plus group here
              Tag_Avail[recap_yr+1,tag_rel,,n_ages,s,sim] <- Tag_Avail[recap_yr+1,tag_rel,,n_ages,s,sim] +
                Tag_Avail[recap_yr,tag_rel,,n_ages,s,sim]  * exp(-tmp_Z[1,,n_ages,s,1])
            } # end else for in plus group
          } # end s loop
        } # end a loop

        # Apply Baranov's to get predicted recaptures
        Pred_Tag_Recap[recap_yr,tag_rel,,,,sim] <- Tag_Reporting[actual_yr,,sim] * (tmp_F / tmp_Z[1,,,,1]) *
          Tag_Avail[recap_yr,tag_rel,,,,sim] * (1 - exp(-tmp_Z[1,,,,1]))

        # Tag Recapture Likelihoods -----------------------------------------------
        # Simulate observed tag recoveries
        # Poisson tag recovery
        for(r in 1:n_regions) {
          for(a in 1:n_ages) {
            for(s in 1:n_sexes) {
              if(tag_like == 0) Obs_Tag_Recap[recap_yr,tag_rel,r,a,s,sim] <- rpois(1, Pred_Tag_Recap[recap_yr,tag_rel,r,a,s,sim]) # Poisson tag recovery
              if(tag_like == 1) Obs_Tag_Recap[recap_yr,tag_rel,r,a,s,sim] <- rnbinom(1, mu = Pred_Tag_Recap[recap_yr,tag_rel,r,a,s,sim], size = tag_nbiom_dispersion[a,s]) # Negbin tag recovery
            } # end s loop
          } # end a loop
        } # end r loop

        # Multinomial tag recovery (release conditioned)
        if(tag_like == 2) {
          tmp_n_tags_rel <- round(sum(Tag_Fish[tag_rel,,,sim])) # Number of initial tags released
          tmp_recap <- aperm(Pred_Tag_Recap[recap_yr,tag_rel,,,,sim, drop = FALSE] / tmp_n_tags_rel, c(4,5,3,1,2,6)) # get recapture probabilities ordered by ages, sexes, regions (i.e., age 1-30, sex 1, region 1, age 1-30, sex 2, region 1, age 1-30, sex 1, region 2, age 1-30, sex 2, region 2 ... ) recapture, tag release, sim ...
          tmp_probs <- c(tmp_recap, 1 - sum(tmp_recap)) # concatenate recapture and non recapture probabilities
          tmp_sim_recap <- stats::rmultinom(1, tmp_n_tags_rel, tmp_probs) # simulate multinomial draws here
          tmp_sim_recap <- aperm(array(tmp_sim_recap[-length(tmp_sim_recap)], dim(tmp_recap)), c(4,5,3,1,2,6)) # remove last group (not recaptured) and then reshape into correct format
          Obs_Tag_Recap[recap_yr,tag_rel,,,,sim] <- tmp_sim_recap # input recaptures from multinomial into observed array
        } # end if for multinomial likelihood (release conditioned)

        # Multinomial tag recovery (recovery conditioned)
        if(tag_like == 3) { # Tag reporting doesn't matter here (if spatially invariant) since it appears in the denominator so it cancels out (when calculating proportion of recaptures)
          tmp_n_tags_recap <- sum(Pred_Tag_Recap[recap_yr,tag_rel,,,,sim]) # Get number of tags to simulate
          tmp_probs <- aperm(Pred_Tag_Recap[recap_yr,tag_rel,,,,sim, drop = FALSE] / tmp_n_tags_recap, c(4,5,3,1,2,6)) # get recapture probabilities ordered by ages, sexes, regions (i.e., age 1-30, sex 1, region 1, age 1-30, sex 2, region 1, age 1-30, sex 1, region 2, age 1-30, sex 2, region 2 ... ) recapture, tag release, sim ...
          tmp_sim_recap <- stats::rmultinom(1, tmp_n_tags_recap, tmp_probs) # simulate multinomial draws here
          tmp_sim_recap <- aperm(array(tmp_sim_recap, dim(tmp_probs)), c(4,5,3,1,2,6)) # reshape into correct format
          Obs_Tag_Recap[recap_yr,tag_rel,,,,sim] <- tmp_sim_recap # input recaptures from multinomial into observed array
        } # end if for Multinomial likelihood (recovery conditioned)

      } # end recap_yr loop
    } # end tag_rel loop
  } # end sim loop

  # Output simulation outputs as a list
  sim_out <- list(init_F = init_F,
                  Fmort = Fmort,
                  fish_sel = fish_sel,
                  ln_rec_devs = ln_rec_devs,
                  M = M,
                  Z = Z,
                  rec_sexratio = rec_sexratio,
                  r0 = r0,
                  WAA = WAA,
                  Maturity_AA = Maturity_AA,
                  init_sigmaR = init_sigmaR,
                  sigmaR = sigmaR,
                  movement_matrix = movement_matrix,
                  Init_NAA = Init_NAA,
                  NAA = NAA,
                  SSB = SSB,
                  Total_Biom = Total_Biom,
                  True_Catch = True_Catch,
                  Obs_Catch = Obs_Catch,
                  CAA = CAA,
                  Obs_FishAgeComps = Obs_FishAgeComps,
                  Obs_SrvIdx = Obs_SrvIdx,
                  True_SrvIdx = True_SrvIdx,
                  Srv_IAA = Srv_IAA,
                  srv_sel = srv_sel,
                  srv_q = srv_q,
                  Obs_SrvAgeComps = Obs_SrvAgeComps,
                  Tag_Release_Ind = as.matrix(tag_rel_indicator),
                  Tag_Reporting = Tag_Reporting,
                  Tag_Fish = Tag_Fish,
                  Tag_Ind_Mort = Tag_Ind_Mort,
                  Tag_Shed = Tag_Shed,
                  Tag_Avail = Tag_Avail,
                  Pred_Tag_Recap = Pred_Tag_Recap,
                  Obs_Tag_Recap = Obs_Tag_Recap
                  )

  saveRDS(sim_out, file = output_path)

}
