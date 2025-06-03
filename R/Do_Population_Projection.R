#' Title Do Population Projections
#'
#' @param n_proj_yrs Number of projection years
#' @param n_regions Number of regions
#' @param n_ages Number of ages
#' @param n_sexes Number of sexes
#' @param sexratio Recruitment sex ratio
#' @param n_fish_fleets Number of fishery fleets
#' @param do_recruits_move Whether recruits move (0 == don't move, 1 == move)
#' @param recruitment Recruitment matrix dimensioned by n_regions, and n_yrs that we want to summarize across, or condition our projection on
#' @param terminal_NAA Terminal Numbers at Age dimensioned by n_regions, n_ages, n_sexes
#' @param terminal_F Terminal fishing mortality rate, dimensioned by n_regions, n_fish_fleets
#' @param natmort Natural mortality, dimensioned by n_regions, n_proj_yrs, n_ages, n_sexes
#' @param WAA Weight at age, dimensioned by n_regions, n_proj_yrs, n_ages, n_sexes
#' @param MatAA Maturity at age, dimensioned by n_regions, n_proj_yrs, n_ages, n_sexes
#' @param fish_sel Fishery selectivity, dimensioned by n_regions, n_proj_yrs, n_ages, n_sexes, n_fish_fleets
#' @param Movement Movement, dimensioned by n_regions, n_regions, n_proj_yrs, n_ages, n_sexes
#' @param f_ref_pt Fishing mortality reference point dimensioned by n_regions and n_proj_yrs
#' @param b_ref_pt Biological reference point dimensioned by n_regions and n_proj_yrs
#' @param HCR_function Function describing a harvest control rule. The function should always have the following arguments: x, which represents SSB, frp, which takes inputs of fishery reference points, and brp, which takes inputs of biological reference points. Any additional arguments should be specified with defaults or hard coded / fixed within the function. e
#' @param recruitment_opt Recruitment simulation option, where options are "inv_gauss", which simulates future recruitment based on the the recruitment values supplied using an inverse gaussian distribution, "mean_rec", which takes the mean of the recruitment values supplied for a given region, and "zero", which assumes that future recruitment does not occur
#' @param fmort_opt Fishing Mortality option, which includes "HCR", which modifies the F reference point using a user supplied HCR_function, or "Input", which uses projected F values supplied by the user.
#' @param t_spawn Fraction time of spawning used to compute projected SSB
#'
#' @returns A list containing projected F, catch, SSB, and Numbers at Age. (Objects are generally dimensioned in the following order: n_regions, n_yrs, n_ages, n_sexes, n_fleets)
#' @export Do_Population_Projection
#'
#' @examples
#' \dontrun{
#' rep <- obj$report(obj$env$last.par.best) # need to have an RTMB object first
#' # Setup necessary inputs
#' n_sims <- 1000
#' t_spawn <- 0
#' sexratio <- 0.5
#' n_proj_yrs <- 15
#' n_regions <- 1
#' n_ages <- length(data$ages)
#' n_sexes <- data$n_sexes
#' n_fish_fleets <- 2
#' do_recruits_move <- 0
#' terminal_NAA <- array(obj$rep$NAA[,length(data$years),,], dim = c(n_regions, n_ages, n_sexes))
#' WAA <- array(rep(data$WAA[,length(data$years),,], each = n_proj_yrs), dim = c(n_regions, n_proj_yrs, n_ages, n_sexes)) # weight at age
#' MatAA <- array(rep(data$MatAA[,length(data$years),,], each = n_proj_yrs), dim = c(n_regions, n_proj_yrs, n_ages, n_sexes)) # maturity at age
#' fish_sel <- array(rep(obj$rep$fish_sel[,length(data$years),,,], each = n_proj_yrs), dim = c(n_regions, n_proj_yrs, n_ages, n_sexes, n_fish_fleets)) # selectivity
#' Movement <- array(rep(obj$rep$Movement[,,length(data$years),,], each = n_proj_yrs), dim = c(n_regions, n_regions, n_proj_yrs, n_ages, n_sexes))
#' terminal_F <- array(obj$rep$Fmort[,length(data$years),], dim = c(n_regions, n_fish_fleets))
#' natmort <- array(obj$rep$natmort[,length(data$years),,], dim = c(n_regions, n_proj_yrs, n_ages, n_sexes))
#' recruitment <- array(obj$rep$Rec[,20:(length(data$years) - 2)], dim = c(n_regions, length(20:length(data$years))))
#'
#' # Define reference points
#' spr_35 <- Get_Reference_Points(data = data,
#'                                rep = rep,
#'                                SPR_x = 0.35, t_spwn = 0, sex_ratio_f = 0.5,
#'                                calc_rec_st_yr = 20, rec_age = 2)
#'
#' spr_40 <- Get_Reference_Points(data = data,
#'                                rep = rep,
#'                                SPR_x = 0.4, t_spwn = 0, sex_ratio_f = 0.5,
#'                                calc_rec_st_yr = 20, rec_age = 2)
#'
#' spr_60 <- Get_Reference_Points(data = data,
#'                                rep = rep,
#'                                SPR_x = 0.6, t_spwn = 0, sex_ratio_f = 0.5,
#'                                calc_rec_st_yr = 20, rec_age = 2)
#'
#' # Extract reference points
#' b40 <- spr_40$b_ref_pt
#' b60 <- spr_60$b_ref_pt
#' f40 <- spr_40$f_ref_pt
#' f35 <- spr_35$f_ref_pt
#' f60 <- spr_60$f_ref_pt
#' # Define the F used for each scenario (Based on BSAI Intro Report)
#' proj_inputs <- list(
#'   # Scenario 1 - Using HCR to adjust maxFABC
#'   list(f_ref_pt = array(f40, dim = c(n_regions, n_proj_yrs)),
#'        b_ref_pt = array(b40, dim = c(n_regions, n_proj_yrs)),
#'        fmort_opt = 'HCR'
#'   ),
#'   # Scenario 2 - Using HCR to adjust maxFABC based on last year's value (constant fraction - author specified F)
#'   list(f_ref_pt = array(f40 * (f40 / 0.086), dim = c(n_regions, n_proj_yrs)),
#'        b_ref_pt = array(b40, dim = c(n_regions, n_proj_yrs)),
#'        fmort_opt = 'HCR'
#'   ),
#'   # Scenario 3 - Using an F input of last 5 years average F, and
#'   list(f_ref_pt = array(mean(rowSums(sabie_rtmb_model$rep$Fmort[1, 60:64, ])), dim = c(n_regions, n_proj_yrs)),
#'        b_ref_pt = NULL,
#'        fmort_opt = 'Input'
#'   ),
#'   # Scenario 4 - Using HCR to adjust F60
#'   list(f_ref_pt = array(f60, dim = c(n_regions, n_proj_yrs)),
#'        b_ref_pt = array(b40, dim = c(n_regions, n_proj_yrs)),
#'        fmort_opt = 'HCR'
#'   ),
#'   # Scenario 5 - F is set at 0
#'   list(f_ref_pt = array(0, dim = c(n_regions, n_proj_yrs)),
#'        b_ref_pt = NULL,
#'        fmort_opt = 'Input'
#'   ),
#'   # Scenario 6 - Using HCR to adjust FOFL
#'   list(f_ref_pt = array(f35, dim = c(n_regions, n_proj_yrs)),
#'        b_ref_pt = array(b40, dim = c(n_regions, n_proj_yrs)),
#'        fmort_opt = 'HCR'
#'   ),
#'   # Scenario 7 - Using HCR to adjust FABC in first 2 projection years, and then later years are adjusting FOFL
#'   list(f_ref_pt = array(c(rep(f40, 2), rep(f35, n_proj_yrs - 2)), dim = c(n_regions, n_proj_yrs)),
#'        b_ref_pt = array(b40, dim = c(n_regions, n_proj_yrs)),
#'        fmort_opt = 'HCR'
#'   )
#' )
#'
#' # store outputs
#' all_scenarios_f <- array(0, dim = c(n_regions, n_proj_yrs, n_sims, length(proj_inputs)))
#' all_scenarios_ssb <- array(0, dim = c(n_regions, n_proj_yrs, n_sims, length(proj_inputs)))
#' all_scenarios_catch <- array(0, dim = c(n_regions, n_proj_yrs, n_fish_fleets, n_sims, length(proj_inputs)))
#'
#' for (i in seq_along(proj_inputs)) {
#'   for (sim in 1:n_sims) {
#'
#'     # do population projection
#'     out <- Do_Population_Projection(n_proj_yrs = n_proj_yrs,
#'                                     n_regions = n_regions,
#'                                     n_ages = n_ages,
#'                                     n_sexes = n_sexes,
#'                                     sexratio = sexratio,
#'                                     n_fish_fleets = n_fish_fleets,
#'                                     do_recruits_move = do_recruits_move,
#'                                     recruitment = recruitment,
#'                                     terminal_NAA = terminal_NAA,
#'                                     terminal_F = terminal_F,
#'                                     natmort = natmort,
#'                                     WAA = WAA,
#'                                     MatAA = MatAA,
#'                                     fish_sel = fish_sel,
#'                                     Movement = Movement,
#'                                     f_ref_pt = proj_inputs[[i]]$f_ref_pt,
#'                                     b_ref_pt = proj_inputs[[i]]$b_ref_pt,
#'                                     HCR_function = HCR_function,
#'                                     recruitment_opt = "inv_gauss",
#'                                     fmort_opt = proj_inputs[[i]]$fmort_opt,
#'                                     t_spawn = t_spawn
#'     )
#'
#'     all_scenarios_ssb[,,sim,i] <- out$proj_SSB
#'     all_scenarios_catch[,,,sim,i] <- out$proj_Catch
#'     all_scenarios_f[,,sim,i] <- out$proj_F[,-(n_proj_yrs+1)] # remove last year, since it's not used
#'   } # end sim loop
#'   print(i)
#' } # end i loop
#' }
Do_Population_Projection <- function(n_proj_yrs = 2,
                                     n_regions,
                                     n_ages,
                                     n_sexes,
                                     sexratio,
                                     n_fish_fleets,
                                     do_recruits_move = 0,
                                     recruitment,
                                     terminal_NAA,
                                     terminal_F,
                                     natmort,
                                     WAA,
                                     MatAA,
                                     fish_sel,
                                     Movement,
                                     f_ref_pt = NULL,
                                     b_ref_pt = NULL,
                                     HCR_function = NULL,
                                     recruitment_opt = "inv_gauss",
                                     fmort_opt = 'HCR',
                                     t_spawn
                                     ) {

# Define Containers -------------------------------------------------------
  fratio <- terminal_F / apply(terminal_F, 1, sum)
  proj_NAA <- array(0, dim = c(n_regions, n_proj_yrs + 1, n_ages, n_sexes))
  proj_ZAA <- array(0, dim = c(n_regions, n_proj_yrs + 1, n_ages, n_sexes))
  proj_FAA <- array(0, dim = c(n_regions, n_proj_yrs + 1, n_ages, n_sexes, n_fish_fleets))
  proj_CAA <- array(0, dim = c(n_regions, n_proj_yrs, n_ages, n_sexes, n_fish_fleets))
  proj_Catch <- array(0, dim = c(n_regions, n_proj_yrs, n_fish_fleets))
  proj_SSB <- array(0, dim = c(n_regions, n_proj_yrs))
  proj_F <- array(0, dim = c(n_regions, n_proj_yrs + 1))

# Start Projection --------------------------------------------------------
  # Input terminal year assessment at age
  proj_NAA[,1,,] <- terminal_NAA

  for(y in 1:n_proj_yrs) {

# Construct Mortality Processes -------------------------------------------

    # use terminal F in the first year (subsequent years use F derived from reference points and HCR)
    if(y == 1) proj_F[,y] <- sum(terminal_F)

    for(f in 1:n_fish_fleets) {
      for(a in 1:n_ages) {
        for(s in 1:n_sexes) {
          # get fishing mortality at age
          proj_FAA[,y,a,s,f] <- proj_F[,y] * fratio[,f] * fish_sel[,y,a,s,f]
        } # end s loop
      } # end a loop
    } # end f loop

    # Get Total Mortality at Age
    proj_ZAA[,y,,] <- natmort[,y,,] + apply(proj_FAA[,y,,,,drop = FALSE], c(1:4), sum)[,1,,] # M and sum F across fleets

# Recruitment Processes ---------------------------------------------------
    tmp_rec <- vector() # temporary variable

    # Use recruits estimated from terminal year in first projection year
    if(y == 1) tmp_rec <- terminal_NAA[,1,]

    if(y > 1) {

      # Inverse Gaussian Recruitment
      if(recruitment_opt == 'inv_gauss') {
        for(r in 1:n_regions) {
          tmp_rec[r] <- rinvgauss_rec(1, recruitment[r,]) # generate inverse gaussian draws
          proj_NAA[r,y,1,] <- tmp_rec[r] * sexratio # input into projected NAA
        } # end r loop
      } # end if

      # Mean Recruitment
      if(recruitment_opt == "mean_rec") {
        for(r in 1:n_regions) {
          tmp_rec[r] <- mean(recruitment[r,]) # get mean recruitment
          proj_NAA[r,y,1,] <- tmp_rec[r] * sexratio # input into projected NAA
        } # end r loop
      }

      # Zero Recruitment
      if(recruitment_opt == "zero") {
        for(r in 1:n_regions) {
          proj_NAA[r,y,1,] <- 0 # input into projected NAA
        } # end r loop
      }

    }

# Movement Processes ------------------------------------------------------
    if(n_regions > 1) {
      # Recruits don't move
      if(do_recruits_move == 0) {
        # Apply movement after ageing processes - start movement at age 2
        for(a in 2:n_ages) for(s in 1:n_sexes) proj_NAA[,y,a,s] = t(proj_NAA[,y,a,s]) %*% Movement[,,y,a,s]
        for(r in 1:n_regions) proj_NAA[r,y,1,] = tmp_rec[r] * sexratio
      } # end if recruits don't move
      # Recruits move here
      if(do_recruits_move == 1) for(a in 1:n_ages) for(s in 1:n_sexes) proj_NAA[,y,a,s] = t(proj_NAA[,y,a,s]) %*% Movement[,,y,a,s]
    } # only compute if spatial

# Mortality and Ageing ----------------------------------------------------
    proj_NAA[,y+1,2:n_ages,] = proj_NAA[,y,1:(n_ages-1),] * exp(-proj_ZAA[,y,1:(n_ages-1),]) # Exponential mortality for individuals not in plus group
    proj_NAA[,y+1,n_ages,] = proj_NAA[,y+1,n_ages,] + proj_NAA[,y,n_ages,] * exp(-proj_ZAA[,y,n_ages,]) # Accumulate plus group

# Derive Biomass ----------------------------------------------------------
    proj_SSB[,y] = apply(proj_NAA[,y,,1,drop = FALSE] * exp(-proj_ZAA[,y,,1,drop = FALSE] * t_spawn) * WAA[,y,,1,drop = FALSE] * MatAA[,y,,1,drop = FALSE], 1, sum) # Spawning Stock Biomass
    if(n_sexes == 1) proj_SSB[,y] = proj_SSB[,y] * 0.5 # If single sex model, multiply SSB calculations by 0.5

# Derive Catches ----------------------------------------------------------
    for(r in 1:n_regions) {
      for(f in 1:n_fish_fleets) {
        for(a in 1:n_ages) {
          for(s in 1:n_sexes) {
            # Get catch at age with Baranov's
            proj_CAA[r,y,a,s,f] <- (proj_FAA[r,y,a,s,f] / proj_ZAA[r,y,a,s]) * proj_NAA[r,y,a,s] * (1 - exp(-proj_ZAA[r,y,a,s]))
          } # end s loop
        } # end a loop

        # Get total catch
        proj_Catch[r,y,f] <- sum(proj_CAA[r,y,,,f] * WAA[r,y,,])

# Project F using HCR and reference points -----------------------------------------------------
        if(fmort_opt == 'HCR') {
          proj_F[r,y+1] <- HCR_function(x = proj_SSB[r,y],
                                        frp = f_ref_pt[r,y],
                                        brp = b_ref_pt[r,y])
        }

# Project F using User Inputs ---------------------------------------------
        if(fmort_opt == 'Input') proj_F[r,y+1] <- f_ref_pt[r,y]

      } # end f loop
    } # end r loop

  } # end y loop

  return(list(proj_F = proj_F,
              proj_Catch = proj_Catch,
              proj_SSB = proj_SSB,
              proj_NAA = proj_NAA,
              proj_ZAA = proj_ZAA)
  )

} # end function
