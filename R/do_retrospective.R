#' Helper function to truncate data years, parameters, and mapping to conduct retrospective diagnostics. Called within do_retrospective function.
#'
#' @param j The years to truncate from the terminal year
#' @param data Data list used for the RTMB model
#' @param parameters Parameter list used for the RTMB model
#' @param mapping Mapping list used for the RTMB model
#'
#' @returns List of data, parameters, and mapping that have truncated dimensions from the original data, parameters, and mapping list
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' retro_list <- retro_truncate_year(j = 0, data, parameters, mapping) # does not remove any data
#' retro_list <- retro_truncate_year(j = 1, data, parameters, mapping) # removes last year of data
#' }
truncate_yr <- function(j,
                        data,
                        parameters,
                        mapping) {

  # set up retro data, parameters, and mapping
  retro_data <- data
  retro_parameters <- parameters
  retro_mapping <- mapping

  # Years
  retro_data$years <- data$years[1:(length(data$years) - j)] # remove j years from years vector
  if(!is.na(sum(retro_data$bias_year))) data$bias_year[4] <- data$bias_year[4] - j # remove j years from bias correction vector

# Recruitment -------------------------------------------------------------

  # Recruitment devs
  retro_parameters$ln_RecDevs <- parameters$ln_RecDevs[,1:(ncol(parameters$ln_RecDevs) - j), drop = FALSE] # Recruitment deviations
  if(any(names(retro_mapping) == 'ln_RecDevs')) retro_mapping$ln_RecDevs <- factor(array(mapping$ln_RecDevs, dim = dim(parameters$ln_RecDevs))[,1:(ncol(parameters$ln_RecDevs) - j), drop = FALSE]) # modify mapping if we have recruitment map


# Fishery -----------------------------------------------------------------

  # Fishery mortality devs
  retro_parameters$ln_F_devs <- parameters$ln_F_devs[,1:(length(data$years) - j),,drop = FALSE] # modify F dev parameters
  retro_mapping$ln_F_devs <- factor(array(mapping$ln_F_devs, dim = dim(parameters$ln_F_devs))[,1:(length(data$years) - j),,drop = FALSE]) # modify map

  # Fishery selectivity deviations
  retro_parameters$ln_fishsel_devs <- parameters$ln_fishsel_devs[,1:(length(data$years) - j),,,,drop = FALSE] # modify parameter length
  retro_mapping$ln_fishsel_devs <- factor(array(mapping$ln_fishsel_devs, dim = dim(parameters$ln_fishsel_devs))[,1:(length(data$years) - j),,,,drop = FALSE]) # modify map
  retro_data$map_ln_fishsel_devs <- data$map_ln_fishsel_devs[,1:(length(data$years) - j),,,,drop = FALSE]

  # Fishery selectivity and catchability blocks
  retro_data$fish_q_blocks <- data$fish_q_blocks[,1:(length(data$years) - j),, drop = FALSE]
  retro_data$fish_sel_blocks <- data$fish_sel_blocks[,1:(length(data$years) - j),, drop = FALSE]

  # Adjust fishery parameter blocks
  retro_parameters$ln_fish_q <- parameters$ln_fish_q[,1:max(retro_data$fish_q_blocks),,drop = FALSE]
  retro_parameters$ln_fish_fixed_sel_pars <- parameters$ln_fish_fixed_sel_pars[,,1:max(retro_data$fish_sel_blocks),,,drop = FALSE]

  # Adjust fishery mapping
  retro_mapping$ln_fish_q <- factor(array(mapping$ln_fish_q, dim = dim(parameters$ln_fish_q))[,1:max(retro_data$fish_q_blocks),,drop = FALSE])
  retro_mapping$ln_fish_fixed_sel_pars <- factor(array(mapping$ln_fish_fixed_sel_pars, dim = dim(parameters$ln_fish_fixed_sel_pars))[,,1:max(retro_data$fish_sel_blocks),,,drop = FALSE])

# Survey ------------------------------------------------------------------

  # Survey selectivity deviations
  retro_parameters$ln_srvsel_devs <- parameters$ln_srvsel_devs[,1:(length(data$years) - j),,,,drop = FALSE] # Survey selectivity deviations
  retro_mapping$ln_srvsel_devs <- factor(array(mapping$ln_srvsel_devs, dim = dim(parameters$ln_srvsel_devs))[,1:(length(data$years) - j),,,,drop = FALSE]) # modify map
  retro_data$map_ln_srvsel_devs <- data$map_ln_srvsel_devs[,1:(length(data$years) - j),,,,drop = FALSE]

  # Survey selectivity and catchability blocks
  retro_data$srv_q_blocks <- data$srv_q_blocks[,1:(length(data$years) - j),, drop = FALSE]
  retro_data$srv_sel_blocks <- data$srv_sel_blocks[,1:(length(data$years) - j),, drop = FALSE]

  # Adjust survey parameter blocks
  retro_parameters$ln_srv_q <- parameters$ln_srv_q[,1:max(retro_data$srv_q_blocks),,drop = FALSE]
  retro_parameters$ln_srv_fixed_sel_pars <- parameters$ln_srv_fixed_sel_pars[,,1:max(retro_data$srv_sel_blocks),,,drop = FALSE]

  # Adjust survey mapping
  retro_mapping$ln_srv_q <- factor(array(mapping$ln_srv_q, dim = dim(parameters$ln_srv_q))[,1:max(retro_data$srv_q_blocks),,drop = FALSE])
  retro_mapping$ln_srv_fixed_sel_pars <- factor(array(mapping$ln_srv_fixed_sel_pars, dim = dim(parameters$ln_srv_fixed_sel_pars))[,,1:max(retro_data$srv_sel_blocks),,,drop = FALSE])

# Movement ----------------------------------------------------------------

  if(data$n_regions > 1) {
    # Movement stuff
    retro_parameters$move_pars <- parameters$move_pars[,,1:(length(data$years) - j),,,drop = FALSE]
    retro_parameters$logit_move_devs <- parameters$logit_move_devs[,,1:(length(data$years) - j),,drop = FALSE]
    retro_mapping$move_pars <- factor(array(mapping$move_pars, dim = dim(parameters$move_pars))[,,1:(length(data$years) - j),,,drop = FALSE])
    retro_mapping$logit_move_devs <- factor(array(mapping$logit_move_devs, dim = dim(parameters$logit_move_devs))[,,1:(length(data$years) - j),,drop = FALSE])
    retro_data$map_Movement_Pars <- data$map_Movement_Pars[,,1:(length(data$years) - j),,,drop = FALSE]
    retro_data$Fixed_Movement <- data$Fixed_Movement[,,1:(length(data$years) - j),,,drop = FALSE]
  }

# Tagging -----------------------------------------------------------------

  if(data$UseTagging == 1) {
    # Tag reporting
    retro_data$Tag_Reporting_blocks <- data$Tag_Reporting_blocks[,1:(length(data$years) - j), drop = FALSE]
    if(!is.na(sum(data$Tag_Reporting_blocks))) retro_parameters$Tag_Reporting_Pars <- parameters$Tag_Reporting_Pars[,1:max(data$Tag_Reporting_blocks),drop = FALSE]
    retro_mapping$Tag_Reporting_Pars <- factor(array(mapping$Tag_Reporting_Pars, dim = dim(parameters$Tag_Reporting_Pars))[,1:max(data$Tag_Reporting_blocks),drop = FALSE])
    retro_data$map_Tag_Reporting_Pars <- data$map_Tag_Reporting_Pars[,1:max(data$Tag_Reporting_blocks),drop = FALSE]

    # Tag cohort stuff
    Tag_Release_Ind <- as.matrix(data$tag_release_indicator)
    retro_data$tag_release_indicator <- as.matrix(Tag_Release_Ind[which(Tag_Release_Ind[,2] %in% 1:(length(data$years) - j)), ])
    retro_data$n_tag_cohorts <- nrow(retro_data$tag_release_indicator)
    retro_data$Tagged_Fish <- data$Tagged_Fish[1:nrow(retro_data$tag_release_indicator),,,drop = FALSE] # remove data (not necessary, but helps with computational cost if using tagging)
    retro_data$Obs_Tag_Recap <- data$Obs_Tag_Recap[,1:nrow(retro_data$tag_release_indicator),,,,drop = FALSE] # remove data (not necessary, but helps with computational cost)
  }


# Data Weights, Composition Stuff, and use indicators ------------------------------------------------------------

  # Data weights and composition stuff
  retro_data$Wt_FishAgeComps <- data$Wt_FishAgeComps[,1:(length(data$years) - j),,,drop = FALSE]
  retro_data$Wt_SrvAgeComps <- data$Wt_SrvAgeComps[,1:(length(data$years) - j),,,drop = FALSE]
  retro_data$Wt_FishLenComps <- data$Wt_FishLenComps[,1:(length(data$years) - j),,,drop = FALSE]
  retro_data$Wt_SrvLenComps <- data$Wt_SrvLenComps[,1:(length(data$years) - j),,,drop = FALSE]
  retro_data$FishAgeComps_Type <- data$FishAgeComps_Type[1:(length(data$years) - j),,drop = FALSE]
  retro_data$FishLenComps_Type <- data$FishLenComps_Type[1:(length(data$years) - j),,drop = FALSE]
  retro_data$SrvLenComps_Type <- data$SrvLenComps_Type[1:(length(data$years) - j),,drop = FALSE]
  retro_data$SrvAgeComps_Type <- data$SrvAgeComps_Type[1:(length(data$years) - j),,drop = FALSE]

  # data use indicators
  retro_data$UseFishAgeComps <- data$UseFishAgeComps[,1:(length(data$years) - j),,drop = FALSE]
  retro_data$UseFishIdx <- data$UseFishIdx[,1:(length(data$years) - j),,drop = FALSE]
  retro_data$UseCatch <- data$UseCatch[,1:(length(data$years) - j),,drop = FALSE]
  retro_data$UseFishLenComps <- data$UseFishLenComps[,1:(length(data$years) - j),,drop = FALSE]
  retro_data$UseSrvAgeComps <- data$UseSrvAgeComps[,1:(length(data$years) - j),,drop = FALSE]
  retro_data$UseSrvIdx <- data$UseSrvIdx[,1:(length(data$years) - j),,drop = FALSE]
  retro_data$UseSrvLenComps <- data$UseSrvLenComps[,1:(length(data$years) - j),,drop = FALSE]


  return(list(retro_data = retro_data,
              retro_parameters = retro_parameters,
              retro_mapping = retro_mapping))
}



#' Functino to run retrospective analyses
#'
#' @param n_retro Number of retrospective peels to do
#' @param data Data list for RTMB model
#' @param parameters Parameter list for RTMB model
#' @param mapping Mapping list for RTMB model
#' @param random Random effects as a character vector - default is NULL
#' @param do_par Whether to do parrallelization, boolean
#' @param do_francis Whether to do francis reweighitng within a given retrospective peel, boolean
#' @param n_francis_iter Number of francis iterations to do
#' @param n_cores Number of cores to use for parrallelization
#' @importFrom stats nlminb optimHess
#' @returns Dataframe of retrospective estiamtes of SSB and recruitment
#' @export do_retrospective
#'
#' @import RTMB
#' @import dplyr
#' @import future.apply
#' @import future
#' @import progressr
#' @importFrom reshape2 melt
#'
#'
#' @examples
#' \dontrun{
#'  # Do retrospective here
#'  ret <- do_retrospective(n_retro = 7, data, parameters, mapping, random = NULL, do_par = TRUE, n_cores = 7, do_francis = TRUE, n_francis_iter = 5)
#'  ggplot(ret, aes(x = Year + 1959, y = value, group = peel, color = 2024 - peel)) +
#'    geom_line(lwd = 1.3) +
#'    facet_wrap(~Type) +
#'    guides (color = guide_colourbar(barwidth = 10, barheight = 1.3)) +
#'    labs(x = 'Year', y = 'Value', color = 'Retrospective Year') +
#'    scale_color_viridis_c() +
#'    theme_bw(base_size = 15) +
#'    theme(legend.position = 'top')
#'
#'  ret %>%
#'    dplyr::mutate(Year = Year + 1959, terminal = 2024 - peel, cohort = Year - 2, years_est = terminal-Year) %>%
#'    filter(Type == 'Recruitment', cohort %in% c(2014:2022), terminal != Year) %>%
#'    ggplot(aes(x = years_est - 1, y = value, group = Year, color = factor(cohort))) +
#'    geom_line(lwd = 1.3) +
#'    geom_point(size = 4) +
#'    theme_bw(base_size = 15) +
#'    labs(x = 'Years since cohort was last estimated', y = 'Recruitment (millions)', color = 'Cohort')
#' }
do_retrospective <- function(n_retro,
                             data,
                             parameters,
                             mapping,
                             random = NULL,
                             do_par,
                             n_cores,
                             do_francis,
                             n_francis_iter = NULL
                             ) {

  # Loop through retrospective (no parrallelization)
  if(do_par == FALSE) {

    retro_all <- data.frame()

    for(j in 0:n_retro) {

      # truncate data
      init <- truncate_yr(j = j, data = data, parameters = parameters, mapping = mapping)

      if(do_francis == FALSE) { # don't do francis within retrospective loop

        # make AD model function
        SPoCK_rtmb_model <- RTMB::MakeADFun(cmb(SPoCK_rtmb, init$retro_data), parameters = init$retro_parameters, map = init$retro_mapping, random = random, silent = T)

        # Now, optimize the function
        SPoCK_optim <- stats::nlminb(SPoCK_rtmb_model$par, SPoCK_rtmb_model$fn, SPoCK_rtmb_model$gr,
                                     control = list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15))
        # newton steps
        try_improve <- tryCatch(expr =
                                  for(i in 1:3) {
                                    g = as.numeric(SPoCK_rtmb_model$gr(SPoCK_optim$par))
                                    h = optimHess(SPoCK_optim$par, fn = SPoCK_rtmb_model$fn, gr = SPoCK_rtmb_model$gr)
                                    SPoCK_optim$par = SPoCK_optim$par - solve(h,g)
                                    SPoCK_optim$objective = SPoCK_rtmb_model$fn(SPoCK_optim$par)
                                  }
                                , error = function(e){e}, warning = function(w){w})

        rep <- SPoCK_rtmb_model$report(SPoCK_rtmb_model$env$last.par.best) # Get report

      } else {

        for(f in 1:n_francis_iter) {

          if(f == 1) { # reset weights at 1 if at the first iteration
            init$retro_data$Wt_FishAgeComps[] <- 1
            init$retro_data$Wt_FishLenComps[] <- 1
            init$retro_data$Wt_SrvAgeComps[] <- 1
            init$retro_data$Wt_SrvLenComps[] <- 1
          } else {
            # get new weights
            wts <- do_francis_reweighting(data = init$retro_data, rep = rep, age_labels = init$retro_data$ages, len_labels = init$retro_data$lens, year_labels = init$retro_data$years)
            init$retro_data$Wt_FishAgeComps[] <- wts$new_fish_age_wts
            init$retro_data$Wt_FishLenComps[] <- wts$new_fish_len_wts
            init$retro_data$Wt_SrvAgeComps[] <- wts$new_srv_age_wts
            init$retro_data$Wt_SrvLenComps[] <- wts$new_srv_len_wts
          }

          # make AD model function
          SPoCK_rtmb_model <- RTMB::MakeADFun(cmb(SPoCK_rtmb, init$retro_data), parameters = init$retro_parameters, map = init$retro_mapping, random = random, silent = T)

          # Now, optimize the function
          SPoCK_optim <- stats::nlminb(SPoCK_rtmb_model$par, SPoCK_rtmb_model$fn, SPoCK_rtmb_model$gr,
                                       control = list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15))
          # newton steps
          try_improve <- tryCatch(expr =
                                    for(i in 1:3) {
                                      g = as.numeric(SPoCK_rtmb_model$gr(SPoCK_optim$par))
                                      h = optimHess(SPoCK_optim$par, fn = SPoCK_rtmb_model$fn, gr = SPoCK_rtmb_model$gr)
                                      SPoCK_optim$par = SPoCK_optim$par - solve(h,g)
                                      SPoCK_optim$objective = SPoCK_rtmb_model$fn(SPoCK_optim$par)
                                    }
                                  , error = function(e){e}, warning = function(w){w})

          rep <- SPoCK_rtmb_model$report(SPoCK_rtmb_model$env$last.par.best) # Get report

        } # end f for francis iteration
      } # end else

      # get ssb and recruitment
      retro_tmp <- reshape2::melt(rep$SSB) %>% dplyr::rename(Region = Var1, Year = Var2) %>% dplyr::mutate(Type = "SSB") %>%
        bind_rows(reshape2::melt(rep$Rec) %>% dplyr::rename(Region = Var1, Year = Var2) %>% dplyr::mutate(Type = "Recruitment")) %>% dplyr::mutate(peel = j)

      retro_all <- rbind(retro_all, retro_tmp) # bind all rows

    } # end j
  } # iterative loop


  # Parrallelize Retrospective Loop
  if(do_par == TRUE) {

    future::plan(future::multisession, workers = n_cores) # set up cores

    progressr::with_progress({

      p <- progressr::progressor(along = 0:n_retro) # progress bar

      retro_all <- future.apply::future_lapply(0:n_retro, function(j) {

        init <- truncate_yr(j = j, data = data, parameters = parameters, mapping = mapping)

        if(do_francis == FALSE) { # don't do francis within retrospective loop

          # make AD model function
          SPoCK_rtmb_model <- fit_model(init$retro_data,
                                        init$retro_parameters,
                                        init$retro_mapping,
                                        random = random,
                                        newton_loops = 3,
                                        silent = T
                                        )

          rep <- SPoCK_rtmb_model$report(SPoCK_rtmb_model$env$last.par.best) # Get report

        } else {

          for(f in 1:n_francis_iter) {

            if(f == 1) { # reset weights at 1 if at the first iteration
              init$retro_data$Wt_FishAgeComps[] <- 1
              init$retro_data$Wt_FishLenComps[] <- 1
              init$retro_data$Wt_SrvAgeComps[] <- 1
              init$retro_data$Wt_SrvLenComps[] <- 1
            } else {
              # get new weights
              wts <- do_francis_reweighting(data = init$retro_data, rep = rep, age_labels = init$retro_data$ages, len_labels = init$retro_data$lens, year_labels = init$retro_data$years)
              init$retro_data$Wt_FishAgeComps[] <- wts$new_fish_age_wts
              init$retro_data$Wt_FishLenComps[] <- wts$new_fish_len_wts
              init$retro_data$Wt_SrvAgeComps[] <- wts$new_srv_age_wts
              init$retro_data$Wt_SrvLenComps[] <- wts$new_srv_len_wts
            }

            # make AD model function
            SPoCK_rtmb_model <- fit_model(init$retro_data,
                                          init$retro_parameters,
                                          init$retro_mapping,
                                          random = random,
                                          newton_loops = 3,
                                          silent = T
            )

            rep <- SPoCK_rtmb_model$report(SPoCK_rtmb_model$env$last.par.best) # Get report

          } # end f for francis iteration
        } # end else

        retro_tmp <- reshape2::melt(rep$SSB) %>% dplyr::rename(Region = Var1, Year = Var2) %>% dplyr::mutate(Type = "SSB") %>%
          bind_rows(reshape2::melt(rep$Rec) %>% dplyr::rename(Region = Var1, Year = Var2) %>% dplyr::mutate(Type = "Recruitment")) %>% dplyr::mutate(peel = j)

        p() # update progress

        retro_tmp

      }, future.seed = TRUE) %>% bind_rows() # bine rows to combine results

      future::plan(future::sequential)  # Reset

    })
  } # do parrallelization for retrospective loop

  return(retro_all)
} # end function

#' Derive relative difference from terminal year from a retrospective analysis.
#'
#' @param retro_data Dataframe outputted from do_retrospective function
#'
#' @returns Returns a data frame with relative difference of SSB and recruitment from the terminal year
#' @export get_retrospective_relative_difference
#'
#' @import dplyr
#' @importFrom tidyr pivot_longer pivot_wider
#' @examples
#' \dontrun{
#'  ret <- do_retrospective(7, data, parameters, mapping, random = NULL, do_par = TRUE, n_cores = 7, do_francis = TRUE, n_francis_iter = 5)
#'  ret_df <- get_retrospective_relative_difference(ret)
#'  ggplot(ret_df %>% filter(Type == 'SSB'), aes(x = Year, y = rd, group = 2024- as.numeric(peel), color = 2024 - as.numeric(peel))) +
#'  geom_hline(yintercept = 0, lty = 2, lwd = 1.3) +
#'    coord_cartesian(ylim = c(-0.4, 0.4)) +
#'    geom_line(lwd = 1.5) +
#'    guides (color = guide_colourbar(barwidth = 15, barheight = 1.3)) +
#'    labs(x = 'Year', y = 'Relative Difference from Terminal Year', color = 'Retrospective Year') +
#'    scale_color_viridis_c() +
#'    theme_bw(base_size = 15) +
#'    theme(legend.position = 'top')
#' }
get_retrospective_relative_difference <- function(retro_data) {

  unique_peels <- length(unique(retro_data$peel)) - 1 # get unique peels

  # Get the terminal year assessment
  terminal <- retro_data %>% dplyr::filter(peel == 0)

  # Get peels
  peels <- retro_data %>% filter(peel != 0) %>%
    tidyr::pivot_wider(names_from = peel, values_from = value)

  # Summarize relative difference
  allret <- terminal %>%
    dplyr::left_join(peels, by = c("Region", "Year", "Type")) %>%
    dplyr::mutate(across(as.character(1:unique_peels), ~ (.x - value) / value, .names = "{.col}"))

  # Pivot longer
  allret <- allret %>%
    dplyr::select(Region, Year, Type, as.character(1:unique_peels)) %>%
    tidyr::pivot_longer(cols = as.character(1:unique_peels), names_to = "peel", values_to = "rd") %>%
    dplyr::mutate(Region = paste("Region", Region))

  return(allret)
}
