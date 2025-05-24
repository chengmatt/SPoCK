#' Title Constructs simulation objects in a new simulation environment for use in simulation functions
#'
#' @param sim_list Simulation list objects
#'
#' @returns A new simulation environment with objects from sim_list
#' @export Setup_sim_env
#'
#' @examples
#' \dontrun{
#' sim_env <- Setup_sim_env(sim_list)
#' }
Setup_sim_env <- function(sim_list) {

  sim_env <- new.env(parent = parent.frame()) # define new environment for simulation

  # Get SPoCK functions in simulation environment
  sim_env$Get_Det_Recruitment <- Get_Det_Recruitment
  sim_env$Get_Tagging_Mortality <- Get_Tagging_Mortality
  sim_env$rdirM <- rdirM

  # output into simulation environment
  list2env(sim_list, envir = sim_env)

  return(sim_env)
}


#' Title Ruin Annual Cycle in Simulation Environment
#'
#' @param y Year index
#' @param sim Simulation index
#' @param sim_env Simulation environment will all the necessary elements to run the annual cycle
#' @export run_annual_cycle
#' @importFrom stats rnorm rmultinom
run_annual_cycle <- function(y,
                             sim,
                             sim_env) {

  # Assign y and sim into simulation environment
  sim_env$y <- y
  sim_env$sim <- sim

  with(sim_env, {
    # Initialize Age Structure ------------------------------------------------
    if(y == 1) {

      # Set up initial equilibrium age structure
      for(r in 1:n_regions) {
        for(s in 1:n_sexes) {
          tmp_cumsum_Z = cumsum(M[r,1,1:(n_ages-1),s,sim] + init_F[r,1,1,sim] * fish_sel[r,1,1:(n_ages-1),s,1,sim]) # cumulative sum of total mortality
          Init_NAA[r,,s,sim] = c(r0[r,1,sim], r0[r,1,sim] * exp(-tmp_cumsum_Z)) * rec_sexratio[r,1,s,sim] # exponential mortality model
        } # end s loop
      } # end r loop

      # Apply annual cycle and iterate to equilibrium
      for(i in 1:init_iter) {
        for(s in 1:n_sexes) {

          Init_NAA_next_year[,1,s,sim] = r0[,1,sim] * rec_sexratio[,1,s,sim] # recruitment

          # recruits don't move
          if(do_recruits_move == 0) for(a in 2:n_ages) Init_NAA[,a,s,sim] = t(Init_NAA[,a,s,sim]) %*% movement_matrix[,,1,a,s,sim] # movement
          # recruits move
          if(do_recruits_move == 1) for(a in 1:n_ages) Init_NAA[,a,s,sim] = t(Init_NAA[,a,s,sim]) %*% movement_matrix[,,1,a,s,sim] # movement

          # ageing and mortality
          Init_NAA_next_year[,2:n_ages,s,sim] = Init_NAA[,1:(n_ages-1),s,sim] * exp(-(M[,1,1:(n_ages-1),s,sim] + (init_F[,1,1,sim] * fish_sel[,1,1:(n_ages-1),s,1,sim])))

          # accumulate plus group
          Init_NAA_next_year[,n_ages,s,sim] = (Init_NAA_next_year[,n_ages,s,sim] * exp(-(M[,1,n_ages,s,sim] + (init_F[,1,1,sim] * fish_sel[,1,n_ages,s,1,sim])))) +
            (Init_NAA[,n_ages,s,sim] * exp(-(M[,1,n_ages,s,sim] + (init_F[,1,1,sim] * fish_sel[,1,n_ages,s,1,sim]))))

          # iterate to next cycle
          Init_NAA = Init_NAA_next_year
        } # end s loop
      } # end i loop

      # Set up initial age deviations
      tmp_ln_init_devs <- NULL # Initialize container vector to allow for global recruitment
      for(r in 1:n_regions) {
        # simulate initial deviations (global recruitment deviations)
        if(init_dd == 0 && is.null(tmp_ln_init_devs)) {
          tmp_ln_init_devs <- stats::rnorm(n_ages-2, 0, init_sigmaR[r,1])
        }
        # simulate initial deviations (local density dependence)
        if(init_dd == 1) {
          tmp_ln_init_devs <- stats::rnorm(n_ages-2, 0, init_sigmaR[r,1])
        }
        # apply deviations
        Init_NAA[r,2:(n_ages-1),s,sim] <- Init_NAA[r,2:(n_ages-1),s,sim] * rep(exp(tmp_ln_init_devs - init_sigmaR[r,1]^2/2), n_sexes)

        # Plug in initial age structure into 1st year (w/o recruitment)
        NAA[r,1,2:n_ages,,sim] <- Init_NAA[r,2:n_ages,s,sim]
      } # end r loop

    } # end initializing age structure

    # Run Annual Cycle --------------------------------------------------------

    ### Recruitment (Year 1 only) ----------------------------------------------
    if(y == 1) {

      tmp_ln_rec_devs <- NULL # Initialize container vector to allow for global recruitment deviations (remains NULL within a given year)

      for(r in 1:n_regions) {

        # Global Recruitment Deviations
        if(rec_dd == 0 && is.null(tmp_ln_rec_devs)) {
          tmp_ln_rec_devs <- stats::rnorm(1, 0, sigmaR[r,y])
        } # Get recruitment deviates (global density dependence)

        # Local Recruitment Deviations
        if(rec_dd == 1) {
          tmp_ln_rec_devs <- ln_rec_devs[r,y,sim] <- stats::rnorm(1, 0, sigmaR[r,y])
        }

        ln_rec_devs[r,y,sim] <- tmp_ln_rec_devs # Input recruitment deviations into vector

        # Get deterministic recruitment
        tmp_det_rec <- Get_Det_Recruitment(recruitment_model = recruitment_opt,
                                           recruitment_dd = rec_dd,
                                           y = y,
                                           rec_lag = rec_lag,
                                           R0 = sum(r0[,y,sim]), # sum to get global R0
                                           Rec_Prop = r0[,y,sim] / sum(r0[,y,sim]), # get R0 proportion
                                           h = h[,y,sim],
                                           n_regions = n_regions,
                                           n_ages = n_ages,
                                           WAA = WAA[,y,,1,sim],
                                           MatAA = Maturity_AA[,y,,1,sim],
                                           natmort = M[,y,,1,sim],
                                           SSB_vals = SSB[,,sim]
        )


        # Input Recruitment into NAA
        for(s in 1:n_sexes) NAA[r,1,1,s,sim] <- tmp_det_rec[r] * exp(ln_rec_devs[r,y,sim] - sigmaR[r,y]^2/2) * rec_sexratio[r,y,s,sim]

      } # end r loop
    } # end if first year recruitment

    for(r in 1:n_regions) {
      for(s in 1:n_sexes) {

        #### Movement ----------------------------------------------------------------
        if(do_recruits_move == 0) for(a in 2:n_ages) for(s in 1:n_sexes) NAA[,y,a,s,sim] = t(NAA[,y,a,s,sim]) %*% movement_matrix[,,y,a,s,sim] # Recruits don't move
        if(do_recruits_move == 1) for(a in 1:n_ages) for(s in 1:n_sexes) NAA[,y,a,s,sim] = t(NAA[,y,a,s,sim]) %*% movement_matrix[,,y,a,s,sim] # Recruits move here

        ### Mortality and Ageing ----------------------------------------------------
        for(a in 1:n_ages) {
          Z[r,y,a,s,sim] <- M[r,y,a,s,sim] + sum(Fmort[r,y,,sim] * fish_sel[r,y,a,s,,sim]) # Total Mortality
          if(a < n_ages) {
            # Exponential mortality for individuals not in plus group (recruits experience mortality )
            NAA[r,y+1,a+1,s,sim] <- NAA[r,y,a,s,sim] * exp(-Z[r,y,a,s,sim])
          } else {
            # Accumulate individuals recently "recruited" into plus group and individuals from previous year
            NAA[r,y+1,n_ages,s,sim] <- NAA[r,y+1,n_ages,s,sim] + NAA[r,y,n_ages,s,sim] * exp(-Z[r,y,a,s,sim])
          } # end else (calculations for plus group)
        } # end a loop
      } # end s loop

      ### Compute Biomass Quantities ----------------------------------------------
      Total_Biom[r,y,sim] <- sum(as.vector(NAA[r,y,,,sim]) * as.vector(WAA[r,y,,,sim])) # Total Biomass
      SSB[r,y,sim] <- sum(as.vector(NAA[r,y,,1,sim]) * as.vector(WAA[r,y,,1,sim]) * Maturity_AA[r,y,,1,sim]) # Spawning Stock Biomass
      if(n_sexes == 1) SSB[r,y,sim] <- SSB[r,y,sim] * 0.5 # If single sex model, multiply SSB calculations by 0.5

      ### Generate Fishery Observations ----------------------------------------------------------
      for(f in 1:n_fish_fleets) {

        # Baranov's catch equation
        CAA[r,y,,,f,sim] <- (Fmort[r,y,f,sim] * fish_sel[r,y,,,f,sim]) / Z[r,y,,,sim] *  NAA[r,y,,,sim] * (1 - exp(-Z[r,y,,,sim]))

        # Generate Catch Data
        True_Catch[r,y,f,sim] <- sum(CAA[r,y,,,f,sim] * WAA[r,y,,,sim]) # True Catch
        Obs_Catch[r,y,f,sim] <- True_Catch[r,y,f,sim] * exp(stats::rnorm(1, 0, sigmaC)) # Observed Catch w/ lognormal deviations

        # Generate Compositions
        for(s in 1:n_sexes) {

          # Split by Region and Sex
          if(comp_strc == 0) {
            tmp_FishAgeComps_Prob <- CAA[r,y,,s,f,sim] # Get CAA vector

            # Simulate multinomial
            if(comp_fish_like[f] == 0) {
              Obs_FishAgeComps[r,y,,s,f,sim] <- array(
                stats::rmultinom(1, ISS_FishAge[r,y,f,sim], prob = tmp_FishAgeComps_Prob / sum(tmp_FishAgeComps_Prob)),
                dim = c(dim(CAA[r,y,,s,f,sim])))
            }

            # Simulate Dirichlet-multinomial
            if(comp_fish_like[f] == 1) {
              Obs_FishAgeComps[r,y,,s,f,sim] <- array(
                rdirM(1, ISS_FishAge[r,y,f,sim], (FishAgeTheta * ISS_FishAge[r,y,f,sim]) * tmp_FishAgeComps_Prob / sum(tmp_FishAgeComps_Prob)),
                dim = dim(tmp_FishAgeComps_Prob))
            }

          } # end if for split approach
        } # end s loop

        # Split by Region, Joint by Sex
        if(comp_strc == 1) {
          tmp_FishAgeComps_Prob <- CAA[r,y,,,f,sim] # Get CAA vector (age and sex vector)

          # Simulate multinomial
          if(comp_fish_like[f] == 0) {
            Obs_FishAgeComps[r,y,,,f,sim] <- array(
              stats::rmultinom(1, ISS_FishAge[r,y,f,sim], prob = tmp_FishAgeComps_Prob / sum(tmp_FishAgeComps_Prob)),
              dim = c(dim(CAA[r,y,,,f,sim, drop = FALSE])))
          }

          # Simulate Dirichlet-multinomial
          if(comp_fish_like[f] == 1) {
            Obs_FishAgeComps[r,y,,,f,sim] <- array(
              rdirM(1, ISS_FishAge[r,y,f,sim], (FishAgeTheta *  ISS_FishAge[r,y,f,sim]) * tmp_FishAgeComps_Prob / sum(tmp_FishAgeComps_Prob)),
              dim = c(dim(CAA[r,y,,,f,sim, drop = FALSE])))
          }

        } # end if joint approach for sexes
      } # end f loop

      ### Generate Survey Observations ----------------------------------------------------------
      for(sf in 1:n_srv_fleets) {

        # Survey Ages Indexed (midpoint year)
        Srv_IAA[r,y,,,sf,sim] <- NAA[r,y,,,sim] * srv_sel[r,y,,,sf,sim] * exp(-0.5 * Z[r,y,,,sim])

        # Generate Survey Index
        True_SrvIdx[r,y,sf,sim] <- srv_q[r,y,sf,sim] * sum(Srv_IAA[r,y,,,sf,sim]) # True Survey Index
        Obs_SrvIdx[r,y,sf,sim] <- srv_q[r,y,sf,sim] * True_SrvIdx[r,y,sf,sim] * exp(stats::rnorm(1, 0, sigmaSrvIdx[r,sf])) # Observed survey index w/ lognormal deviations

        # Generate Survey Compositions
        for(s in 1:n_sexes) {
          # Split Approach
          if(comp_strc == 0) {
            tmp_SrvAgeComps_Prob <- Srv_IAA[r,y,,s,sf,sim] # Get survey index at age vector

            # Simulate Multinomial
            if(comp_srv_like[sf] == 0) {
              Obs_SrvAgeComps[r,y,,s,sf,sim] <- array(
                stats::rmultinom(1, ISS_SrvAge, tmp_SrvAgeComps_Prob / sum(tmp_SrvAgeComps_Prob)),
                dim = c(dim(Srv_IAA[r,y,,s,sf,sim])))
            }

            # Simulate Dirichlet-Multinomial
            if(comp_srv_like[sf] == 1) {
              Obs_SrvAgeComps[r,y,,s,sf,sim] <- array(
                rdirM(1, ISS_SrvAge, (SrvAgeTheta * ISS_SrvAge) * tmp_SrvAgeComps_Prob / sum(tmp_SrvAgeComps_Prob)),
                dim = c(dim(Srv_IAA[r,y,,s,sf,sim])))
            }

          } # end if Split
        } # end s loop

        # Joint Approach
        if(comp_strc == 1) {
          tmp_SrvAgeComps_Prob <- Srv_IAA[r,y,,,sf,sim] # Get survey vector (age and sex vector)

          # Simulate Multinomial
          if(comp_srv_like[sf] == 0) {
            Obs_SrvAgeComps[r,y,,,sf,sim] <- array(
              stats::rmultinom(1, ISS_SrvAge, tmp_SrvAgeComps_Prob / sum(tmp_SrvAgeComps_Prob)),
              dim = c(dim(Srv_IAA[r,y,,,sf,sim,drop=FALSE])))
          }

          # Simulate Dirichlet-Multinomial
          if(comp_srv_like[sf] == 1) {
            Obs_SrvAgeComps[r,y,,,sf,sim] <- array(
              rdirM(1, ISS_SrvAge, (SrvAgeTheta * ISS_SrvAge) * tmp_SrvAgeComps_Prob / sum(tmp_SrvAgeComps_Prob)),
              dim = c(dim(Srv_IAA[r,y,,,sf,sim,drop=FALSE])))
          }

        } # end if Joint by sex
      } # end sf loop


      ### Release Tags ------------------------------------------------
      # Get indices for tag cohorts in the current year and region
      tag_rel <- which(tag_rel_indicator$regions == r & tag_rel_indicator$tag_yrs == y) # Get tag cohort (release event)

      # Release Tags if any events
      if(length(tag_rel) != 0) {
        # Tag Indexing
        tag_rel_region <- tag_rel_indicator$regions[tag_rel] # tag release region
        tag_rel_yr <- tag_rel_indicator$tag_yrs[tag_rel] # tag release year

        # Release Tagged Fish
        n_tags_rel <- round(Obs_SrvIdx[tag_rel_region,tag_rel_yr,1,sim] / sum(Obs_SrvIdx[,tag_rel_yr,1,sim]) * n_tags) # Number of tags apportioned across regions
        tmp_SrvAgeComps_Prob <- as.vector(Srv_IAA[tag_rel_region, tag_rel_yr,,,1,sim]) # Tagged Fish Proportions at age and sex (from first survey)
        tagged_fish <- round((tmp_SrvAgeComps_Prob / sum(tmp_SrvAgeComps_Prob)) * n_tags_rel) # Tagged fish
        Tag_Fish[tag_rel,,,sim] <- array(tagged_fish, dim = c(n_ages, n_sexes)) # Input into Tag_Fish array and apply initial tag induced mortality
      }

    } # end r loop

    #### Generate Tag Recaptures ------------------------------------------------
    for(tag_rel in 1:n_tag_rel_events) {

      # get indexing
      tag_rel_region <- tag_rel_indicator$regions[tag_rel] # tag release region
      tag_rel_yr <- tag_rel_indicator$tag_yrs[tag_rel] # tag release year

      # Skipping stuff if hasn't occurred yet, or if max liberty
      if(tag_rel_yr > y) next # skip if tagging hasn't occurred
      recap_yr <- y - tag_rel_yr + 1 # get tag liberty
      if(recap_yr > max_liberty) next # skip if max liberty

      # Input tagged fish into available tags for recapture and adjust initial number of tagged fish for tag induced mortality (exponential mortality process)
      if(recap_yr == 1) Tag_Avail[1,tag_rel,tag_rel_region,,,sim] <- Tag_Fish[tag_rel,,,sim] * exp(-Tag_Ind_Mort)

      # Get mortality estimates
      tmp_F <- Get_Tagging_Mortality(tag_selex = 5, # F weighted selectivity
                                     tag_natmort = 3, # Age and sex-specific natural mortality
                                     # Reformat fishing mortality
                                     Fmort = array(Fmort[,,,sim,drop = FALSE], dim = c(n_regions, n_yrs, n_fish_fleets)),
                                     # Reformat natural mortality
                                     natmort = array(M[,,,,sim,drop = FALSE], dim = c(n_regions, n_yrs, n_ages, n_sexes)),
                                     Tag_Shed = Tag_Shed,
                                     # Reformat fishery selectivity
                                     fish_sel = array(fish_sel[,,,,,sim,drop = FALSE], dim = c(n_regions, n_yrs, n_ages, n_sexes, n_fish_fleets)),
                                     n_regions = n_regions,
                                     n_ages = n_ages,
                                     n_sexes = n_sexes,
                                     n_fish_fleets = n_fish_fleets,
                                     y = y,
                                     what = "F"
      )

      # Get total mortality
      tmp_Z <- Get_Tagging_Mortality(tag_selex = 5, # F weighted selectivity
                                     tag_natmort = 3, # Age and sex-specific natural mortality
                                     # Reformat fishing mortality
                                     Fmort = array(Fmort[,,,sim,drop = FALSE], dim = c(n_regions, n_yrs, n_fish_fleets)),
                                     # Reformat natural mortality
                                     natmort = array(M[,,,,sim,drop = FALSE], dim = c(n_regions, n_yrs, n_ages, n_sexes)),
                                     Tag_Shed = Tag_Shed,
                                     # Reformat fishery selectivity
                                     fish_sel = array(fish_sel[,,,,,sim,drop = FALSE], dim = c(n_regions, n_yrs, n_ages, n_sexes, n_fish_fleets)),
                                     n_regions = n_regions,
                                     n_ages = n_ages,
                                     n_sexes = n_sexes,
                                     n_fish_fleets = n_fish_fleets,
                                     y = y,
                                     what = "Z"
      )

      # Discount with tagging time (t_tagging) if it doesn't happen at the start of the year
      if(recap_yr == 1 && t_tagging != 1) tmp_Z <- tmp_Z * t_tagging

      # Move tagged fish around
      if(t_tagging != 1 && recap_yr == 1) { # Movement does not occur if tagging does not happen at start of year
      } else{
        for(a in 1:n_ages) {
          for(s in 1:n_sexes) {
            Tag_Avail[recap_yr,tag_rel,,a,s,sim] <- Tag_Avail[recap_yr,tag_rel,,a,s,sim] %*% movement_matrix[,,y,a,s,sim]
          } # end s loop
        } # end a loop
      }

      # Apply mortality and ageing to tagged fish
      for(a in 1:n_ages) {
        for(s in 1:n_sexes) {
          if(a < n_ages) { # If not in plus group
            Tag_Avail[recap_yr+1,tag_rel,,a+1,s,sim] <- Tag_Avail[recap_yr,tag_rel,,a,s,sim] * exp(-tmp_Z[,1,a,s])
          } else{ # Accumulate plus group here
            Tag_Avail[recap_yr+1,tag_rel,,n_ages,s,sim] <- Tag_Avail[recap_yr+1,tag_rel,,n_ages,s,sim] +
              Tag_Avail[recap_yr,tag_rel,,n_ages,s,sim]  * exp(-tmp_Z[,1,n_ages,s])
          } # end else for in plus group
        } # end s loop
      } # end a loop

      # Apply Baranov's to get predicted recaptures
      Pred_Tag_Recap[recap_yr,tag_rel,,,,sim] <- Tag_Reporting[,y,sim] * (tmp_F[,1,,] / tmp_Z[,1,,]) *
        Tag_Avail[recap_yr,tag_rel,,,,sim] * (1 - exp(-tmp_Z[,1,,]))

      # Simulate Tag Recoveries
      if(tag_like %in% c(0,1)) {
        for(r in 1:n_regions) {
          for(a in 1:n_ages) {
            for(s in 1:n_sexes) {

              # Poisson Tag Recovery
              if(tag_like == 0){
                Obs_Tag_Recap[recap_yr,tag_rel,r,a,s,sim] <- rpois(1, Pred_Tag_Recap[recap_yr,tag_rel,r,a,s,sim])
              }

              # Negative Binomial Tag Recovery
              if(tag_like == 1) {
                Obs_Tag_Recap[recap_yr,tag_rel,r,a,s,sim] <- rnbinom(1, mu = Pred_Tag_Recap[recap_yr,tag_rel,r,a,s,sim], size = Tag_Like_Pars)
              }

            } # end s loop
          } # end a loop
        } # end r loop
      } # end if

      # Multinomial tag recovery (release conditioned)
      if(tag_like == 2) {
        # Get number of initial tags released
        tmp_n_tags_rel <- round(sum(Tag_Fish[tag_rel,,,sim]))
        # get recapture probabilities ordered by ages, sexes, regions
        tmp_recap <- aperm(Pred_Tag_Recap[recap_yr,tag_rel,,,,sim, drop = FALSE] / tmp_n_tags_rel, c(4,5,3,1,2,6))
        # concatenate recapture and non recapture probabilities
        tmp_probs <- c(tmp_recap, 1 - sum(tmp_recap))
        # simulate multinomial draws here
        tmp_sim_recap <- stats::rmultinom(1, tmp_n_tags_rel, tmp_probs)
        # remove last group (not recaptured) and then reshape into correct format
        tmp_sim_recap <- aperm(array(tmp_sim_recap[-length(tmp_sim_recap)], dim(tmp_recap)), c(4,5,3,1,2,6))
        # input recaptures from multinomial into observed array
        Obs_Tag_Recap[recap_yr,tag_rel,,,,sim] <- tmp_sim_recap
      } # end if for multinomial likelihood (release conditioned)

      # Multinomial tag recovery (recovery conditioned)
      if(tag_like == 3) {
        # Get number of tags to simulate
        tmp_n_tags_recap <- sum(Pred_Tag_Recap[recap_yr,tag_rel,,,,sim]) # Get number of tags to simulate
        # get recapture probabilities ordered by ages, sexes, regions
        tmp_probs <- aperm(Pred_Tag_Recap[recap_yr,tag_rel,,,,sim, drop = FALSE] / tmp_n_tags_recap, c(4,5,3,1,2,6))
        # simulate multinomial draws here
        tmp_sim_recap <- stats::rmultinom(1, tmp_n_tags_recap, tmp_probs)
        # reshape into correct format
        tmp_sim_recap <- aperm(array(tmp_sim_recap, dim(tmp_probs)), c(4,5,3,1,2,6))
        # input recaptures from multinomial into observed array
        Obs_Tag_Recap[recap_yr,tag_rel,,,,sim] <- tmp_sim_recap
      } # end if for Multinomial likelihood (recovery conditioned)
    } # end tag_rel loop


    ### Compute Recruitment for Next Year ---------------------------------------
    tmp_ln_rec_devs_next <- NULL # Initialize container vector for next year's recruitment deviations
    if(y < n_yrs) {
      for(r in 1:n_regions) {

        # Get Recruitment Deviations for next year
        # Global Recruitment Deviations
        if(rec_dd == 0 && is.null(tmp_ln_rec_devs_next)) {
          tmp_ln_rec_devs_next <- stats::rnorm(1, 0, sigmaR[r,y+1])
        }

        # Local Recruitment Deviations
        if(rec_dd == 1) {
          tmp_ln_rec_devs_next <- ln_rec_devs[r,y+1,sim] <- stats::rnorm(1, 0, sigmaR[r,y+1])
        }

        ln_rec_devs[r,y+1,sim] = tmp_ln_rec_devs_next # Input recruitment deviations into vector

        # Get Deterministic Recruitment for next year
        tmp_det_rec_next <- Get_Det_Recruitment(recruitment_model = recruitment_opt,
                                                recruitment_dd = rec_dd,
                                                y = y+1,
                                                rec_lag = rec_lag,
                                                R0 = sum(r0[,y,sim]), # sum to get global R0
                                                Rec_Prop = r0[,y,sim] / sum(r0[,y,sim]), # get R0 proportion
                                                h = h[,y,sim],
                                                n_regions = n_regions,
                                                n_ages = n_ages,
                                                WAA = WAA[,y,,1,sim],
                                                MatAA = Maturity_AA[,y,,1,sim],
                                                natmort = M[,y,,1,sim],
                                                SSB_vals = SSB[,,sim]
        )

        # Store next year's recruitment (will be used when next year starts)
        for(s in 1:n_sexes)  NAA[r,y+1,1,s,sim] <- tmp_det_rec_next[r] * exp(ln_rec_devs[r,y+1,sim] - sigmaR[r,y+1]^2/2) * rec_sexratio[r,y+1,s,sim]

      } # end r loop
    } # end if

  }) # end simulation environment

  return(invisible(NULL))

}

#' Simulates a static spatial, sex, and age-structured population (no feedback loop)
#'
#' @param output_path path to output simulation objects
#' @param sim_list Simulation list objects
#'
#' @returns a list object with a bunch of simulated values and outputs
#' @export Simulate_Pop_Static
#' @examples
#' \dontrun{
#'   # Set up model dimensions
#'  sim_list <- Setup_Sim_Dim(n_sims = 100,
#'                            n_yrs = 10,
#'                            n_regions = 2,
#'                            n_ages = 8,
#'                            n_sexes = 1,
#'                            n_fish_fleets = 1,
#'                            n_srv_fleets = 1
#'  )
#'
#'  # set up containers
#'  sim_list <- Setup_Sim_Containers(sim_list)
#'
#'  # Setup fishing mortality
#'  sim_list <- Setup_Sim_FishMort(sim_list = sim_list,
#'                                 sigmaC = 1e-3,
#'                                 init_F = matrix(0, nrow = sim_list$n_regions, ncol = sim_list$n_fish_fleets),
#'                                 Fmort_pattern = matrix(c('two-way', "two-way"), nrow = sim_list$n_regions, ncol = sim_list$n_fish_fleets),
#'                                 Fmort_start = matrix(c(0.01, 0.01), nrow = sim_list$n_regions, ncol = sim_list$n_fish_fleets),
#'                                 Fmort_fct = matrix(c(15, 15), nrow = sim_list$n_regions, ncol = sim_list$n_fish_fleets),
#'                                 proc_error = TRUE,
#'                                 proc_error_sd = 0.1)
#'
#'  # Setup fishery selectivity
#'  sim_list <- Setup_Sim_FishSel(sel_model = matrix(c('logistic', "logistic"), nrow = sim_list$n_regions, ncol = sim_list$n_fish_fleets),
#'                                # a50, k for logistic shared across regions
#'                                fixed_fish_sel_pars = array(c(3,3,1,1), dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets, 2)),
#'                                sim_list = sim_list
#'  )
#'
#'  # Setup survey catchability and selectivity
#'  sim_list <- Setup_Sim_Survey(
#'    sim_list = sim_list,
#'    sigmaSrvIdx = array(0.2, dim = c(sim_list$n_regions, sim_list$n_srv_fleets)), # survey observation error
#'    base_srv_q = array(1, dim = c(sim_list$n_regions, sim_list$n_srv_fleets)), # base survey catchability value
#'    srv_q_pattern = matrix(c('constant', "constant"), nrow = sim_list$n_regions, ncol = sim_list$n_srv_fleets), # catchability pattern
#'    sel_model = matrix(c('logistic', "logistic"), nrow = sim_list$n_regions, ncol = sim_list$n_srv_fleets), # selectivity model
#'    # a50, k, for logistic shared across regions
#'    fixed_srv_sel_pars = array(c(3,3,1,1), dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_srv_fleets, 2))
#'  )
#'
#'  # Setup recruitment stuff
#'  sim_list <- Setup_Sim_Rec(
#'    sim_list = sim_list,
#'    do_recruits_move = "dont_move", # == 0, recruits don't move , == 1 recruits move
#'    base_rec_sexratio = 1, # single sex
#'    rec_sexratio_vary = "constant",
#'    base_r0 = c(50, 50),
#'    r0_vary = "constant",
#'    base_h = c(0.8, 0.8),
#'    init_sigmaR = 0.5,
#'    sigmaR = 0.5,
#'    recruitment_opt = "bh_rec",
#'    rec_dd = "global",
#'    init_dd = "global",
#'    rec_lag = 1
#'  )
#'
#'  # Setup biologicals
#'  sim_list <- Setup_Sim_Biologicals(
#'    sim_list = sim_list,
#'    base_M_value = array(0.5, dim = c(sim_list$n_regions, sim_list$n_ages, sim_list$n_sims)),
#'    M_pattern = "constant",
#'    base_WAA_values = array(rep(5 * (1 - exp(-0.1 * 1:sim_list$n_ages)), each = sim_list$n_regions * sim_list$n_sexes),
#'                            dim = c(sim_list$n_regions, sim_list$n_ages, sim_list$n_sexes)),
#'    WAA_pattern = "constant",
#'    base_Maturity_AA_values = array(rep(1 / (1 + exp(-0.3 * 1:sim_list$n_ages)), each = sim_list$n_regions * sim_list$n_sexes),
#'                                    dim = c(sim_list$n_regions, sim_list$n_ages, sim_list$n_sexes)),
#'    Maturity_AA_pattern = "constant"
#'  )
#'
#'  # Setup tagging stuff
#'  sim_list <- Setup_Sim_Tagging(
#'    sim_list = sim_list,
#'    n_tags = 5000,
#'    max_liberty = 30,
#'    tag_years = seq(1, sim_list$n_yrs, 3),
#'    t_tagging = 0.5,
#'    base_Tag_Reporting = c(0.2, 0.2),
#'    Tag_Reporting_pattern = "constant",
#'    Tag_Ind_Mort = 0,
#'    Tag_Shed = 0
#'  )
#'
#'  # Setup observation processes
#'  sim_list <- Setup_Sim_Observation_Proc(
#'    sim_list = sim_list,
#'    Comp_Structure = "spltR_jntS",
#'    Comp_Srv_Like = "Multinomial",
#'    Comp_Fish_Like = "Multinomial",
#'    ISS_FishAge_Pattern = 'F_pattern',
#'    FishAgeTheta = 3,
#'    SrvAgeTheta = 2,
#'    Srv_Like_Pars = NA,
#'    base_ISS_FishAge = 200,
#'    base_ISS_SrvAge = 200,
#'    Tag_Like = "Poisson",
#'    Tag_Like_Pars = NA
#'  )
#'
#'  # IID Movement Matrix across years and ages
#'  ref <- 1
#'  movement_matrix <- array(0, dim = c(sim_list$n_regions, sim_list$n_regions, sim_list$n_yrs, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims)) # From, To
#'  base <- matrix(0, sim_list$n_regions, sim_list$n_regions)
#'
#'  # Plug in movement process error
#'  for(sim in 1:sim_list$n_sims) {
#'    for(a in 1:sim_list$n_ages) {
#'      for(s in 1:sim_list$n_sexes) {
#'        for(y in 1:sim_list$n_yrs) {
#'          for(r in 1:sim_list$n_regions) {
#'            if(a > 1) pe_err <- rnorm(length(tmp_move[-ref]), 0, 0.4) # logit proces error
#'            else pe_err <- 0
#'            tmp_move <- base[r,]
#'            tmp_move[-ref] <- tmp_move[-ref] + pe_err
#'            movement_matrix[r,,y,a,s,sim] <- exp(tmp_move) / sum(exp(tmp_move))
#'          }
#'        } # end y loop
#'      } # end s loop
#'    } # end a loop
#'  } # end sim loop
#'
#'  sim_list$movement_matrix <- movement_matrix
#'
#'  # Run Simulation ----------------------------------------------------------
#'  Simulate_Pop_Static(sim_list = sim_list, output_path = here("sim_out.RDS"))
#' }
Simulate_Pop_Static <- function(sim_list,
                                output_path) {

  # Setup simulation environment
  sim_env <- Setup_sim_env(sim_list)

  # Start Simulation
  for (sim in 1:sim_env$n_sims) {
    for (y in 1:sim_env$n_yrs) {
      # Run annual cycle here
      run_annual_cycle(y, sim, sim_env)
    } # end y loop
  } # end sim loop

  # Output simulation outputs as a list
  sim_out <- list(init_F = sim_env$init_F,
                  Fmort = sim_env$Fmort,
                  fish_sel = sim_env$fish_sel,
                  ln_rec_devs = sim_env$ln_rec_devs,
                  M = sim_env$M,
                  Z = sim_env$Z,
                  rec_sexratio = sim_env$rec_sexratio,
                  r0 = sim_env$r0,
                  WAA = sim_env$WAA,
                  Maturity_AA = sim_env$Maturity_AA,
                  init_sigmaR = sim_env$init_sigmaR,
                  sigmaR = sim_env$sigmaR,
                  movement_matrix = sim_env$movement_matrix,
                  Init_NAA = sim_env$Init_NAA,
                  NAA = sim_env$NAA,
                  SSB = sim_env$SSB,
                  Total_Biom = sim_env$Total_Biom,
                  True_Catch = sim_env$True_Catch,
                  Obs_Catch = sim_env$Obs_Catch,
                  CAA = sim_env$CAA,
                  Obs_FishAgeComps = sim_env$Obs_FishAgeComps,
                  Obs_SrvIdx = sim_env$Obs_SrvIdx,
                  True_SrvIdx = sim_env$True_SrvIdx,
                  Srv_IAA = sim_env$Srv_IAA,
                  srv_sel = sim_env$srv_sel,
                  srv_q = sim_env$srv_q,
                  Obs_SrvAgeComps = sim_env$Obs_SrvAgeComps,
                  Tag_Release_Ind = as.matrix(sim_env$tag_rel_indicator),
                  Tag_Reporting = sim_env$Tag_Reporting,
                  Tag_Fish = sim_env$Tag_Fish,
                  Tag_Ind_Mort = sim_env$Tag_Ind_Mort,
                  Tag_Shed = sim_env$Tag_Shed,
                  Tag_Avail = sim_env$Tag_Avail,
                  Pred_Tag_Recap = sim_env$Pred_Tag_Recap,
                  Obs_Tag_Recap = sim_env$Obs_Tag_Recap
                  )

  # save RDS file
  saveRDS(sim_out, file = output_path)

} # end function

#' Title Truncates data to feedback simulation loop year
#'
#' @param sim_env Simulation environment defined using Simulate_Pop_Feedback
#' @param sim_list Simulation list that defines the operating model options
#' @param y y index used for annual cycle
#' @param sim sim index used for annual cycle
#' @param skeleton_data Skeleton data list for estimation method (SPoCK) where all dimensions are defined using the total number of simulation years
#' @param skeleton_parameters Skeleton parameters list for estimation method (SPoCK) where all dimensions are defined using the total number of simulation years
#' @param skeleton_mapping Skeleton mapping list for estimation method (SPoCK) where all dimensions are defined using the total number of simulation years
#'
#' @export Get_Feedback_Data
#' @returns Returns data, parameters, and mapping list truncated to match the feedback year
Get_Feedback_Data <- function(sim_env,
                              sim_list,
                              y,
                              sim,
                              skeleton_data,
                              skeleton_parameters,
                              skeleton_mapping
                              ) {

  # Get total years in closed loop feedback
  n_yrs <- length(skeleton_data$years)

  # Get skeleton assessment inputs with uninitialized data and truncate arrays
  feedback <- truncate_yr(j = n_yrs - y, data = skeleton_data, parameters = skeleton_parameters, mapping = skeleton_mapping)

  # Input data from simulation list and environment into feedback data (the following code simply subsets data up to the feedback year)

  # Weight at age
  feedback$retro_data$WAA <- array(sim_list$WAA[,1:y,,,sim, drop = FALSE],
                                   dim = c(feedback$retro_data$n_regions, length(feedback$retro_data$years),
                                           length(feedback$retro_data$ages), feedback$retro_data$n_sexes))
  # Maturity at age
  feedback$retro_data$MatAA <- array(sim_list$Maturity_AA[,1:y,,,sim, drop = FALSE],
                                     dim = c(feedback$retro_data$n_regions, length(feedback$retro_data$years),
                                             length(feedback$retro_data$ages), feedback$retro_data$n_sexes))

  # Catch Stuff
  feedback$retro_data$ObsCatch <- array(sim_env$Obs_Catch[,1:y,,sim, drop = FALSE],
                                        dim = c(feedback$retro_data$n_regions, length(feedback$retro_data$years),
                                                feedback$retro_data$n_fish_fleets))

  # Fishery Age Compositions
  feedback$retro_data$ObsFishAgeComps <- array(sim_env$Obs_FishAgeComps[,1:y,,,,sim,drop = FALSE],
                                               dim = c(feedback$retro_data$n_regions, length(feedback$retro_data$years),
                                                       length(feedback$retro_data$ages), feedback$retro_data$n_sexes, feedback$retro_data$n_fish_fleets))

  # Survey Index
  feedback$retro_data$ObsSrvIdx = array(sim_env$Obs_SrvIdx[,1:y,,sim,drop = FALSE],
                                        dim = c(feedback$retro_data$n_regions, length(feedback$retro_data$years),
                                                feedback$retro_data$n_srv_fleets))

  # Survey Age Compositions
  feedback$retro_data$ObsSrvAgeComps = array(sim_env$Obs_SrvAgeComps[,1:y,,,,sim,drop = FALSE],
                                             dim = c(feedback$retro_data$n_regions, length(feedback$retro_data$years),
                                                     length(feedback$retro_data$ages), feedback$retro_data$n_sexes,
                                                     feedback$retro_data$n_srv_fleets))

  # Tagging
  Tag_Release_Ind <- as.matrix(sim_env$tag_rel_indicator)
  Tag_Release_Ind <- Tag_Release_Ind[which(Tag_Release_Ind[,2] %in% 1:y),]
  feedback$retro_data$tag_release_indicator <- Tag_Release_Ind
  feedback$retro_data$Tagged_Fish <- array(sim_env$Tag_Fish[1:nrow(Tag_Release_Ind),,,sim],
                                           dim = c(dim(sim_env$Tag_Fish[1:nrow(Tag_Release_Ind),,,sim]),
                                                   feedback$retro_data$n_sexes))
  feedback$retro_data$Obs_Tag_Recap <- array(sim_env$Obs_Tag_Recap[,1:nrow(Tag_Release_Ind),,,,sim],
                                            dim = c(dim(sim_env$Obs_Tag_Recap[,1:nrow(Tag_Release_Ind),,,,sim]),
                                                    feedback$retro_data$n_sexes))

  return(feedback)

}


