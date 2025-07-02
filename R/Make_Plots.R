#' Get Time Series Plots
#'
#' @param rep List of n_models of `SPoRC` report lists
#' @param sd_rep List of n_models of `SPoRC` sdreport lists
#' @param model_names Vector of model names
#'
#' @returns Plots of spawning biomass, total biomass, recruitment, and fishing mortality time-series across models
#' @export get_ts_plot
#'
#' @examples
#' \dontrun{
#'   get_ts_plot(list(rep1, rep2), list(sd_rep1, sd_rep2), c("Model1", "Model2"))
#' }
get_ts_plot <- function(rep,
                        sd_rep,
                        model_names
                         ) {

  biom_rec_df <- data.frame() # empty dataframe to bind

  for(i in 1:length(rep)) {

    # Spawning Stock Biomass
    ssb_plot_df <- reshape2::melt(rep[[i]]$SSB) %>%
      dplyr::rename(Region = Var1, Year = Var2) %>%
      dplyr::bind_cols(se = sd_rep[[i]]$sd[names(sd_rep[[i]]$value) == "log(SSB)"]) %>%
      dplyr::mutate(lwr = exp(log(value) - 1.96 * se),
                    upr = exp(log(value) + 1.96 * se),
                    Region = paste("Region", Region),
                    Type = 'SSB',
                    Model = model_names[i])

    # Total Biomass
    totbiom_plot_df <- reshape2::melt(rep[[i]]$Total_Biom) %>%
      dplyr::rename(Region = Var1, Year = Var2) %>%
      dplyr::bind_cols(se = sd_rep[[i]]$sd[names(sd_rep[[i]]$value) == "log(Total_Biom)"]) %>%
      dplyr::mutate(lwr = exp(log(value) - 1.96 * se),
                    upr = exp(log(value) + 1.96 * se),
                    Region = paste("Region", Region),
                    Type = 'Total Biom',
                    Model = model_names[i])

    # Recruitment
    rec_plot_df <- reshape2::melt(rep[[i]]$Rec) %>%
      dplyr::rename(Region = Var1, Year = Var2) %>%
      dplyr::bind_cols(se = sd_rep[[i]]$sd[names(sd_rep[[i]]$value) == "log(Rec)"]) %>%
      dplyr::mutate(lwr = exp(log(value) - 1.96 * se),
                    upr = exp(log(value) + 1.96 * se),
                    Region = paste("Region", Region),
                    Type = 'Recruitment',
                    Model = model_names[i])

    # Fishing Mortality
    f_plot_df <- reshape2::melt(rep[[i]]$Fmort) %>%
      dplyr::rename(Region = Var1, Year = Var2, Type = Var3) %>%
      dplyr::mutate(Region = paste("Region", Region),
                    Type = paste("Fleet", Type, "F"),
                    lwr = NA,
                    upr = NA,
                    se = NA,
                    Model = model_names[i])

    # bind together
    biom_rec_df <- rbind(ssb_plot_df, totbiom_plot_df, rec_plot_df, f_plot_df, biom_rec_df)
  }

  # Plot time series
  ts_plot <- ggplot2::ggplot(biom_rec_df,
                             ggplot2::aes(x = Year, y = value, ymin = lwr, ymax = upr, color = factor(Model), fill = factor(Model))) +
    ggplot2::geom_line(lwd = 0.9) +
    ggplot2::geom_ribbon(alpha = 0.3, color = NA) +
    ggplot2::facet_grid(Type~Region, scales = 'free') +
    ggplot2::labs(x = 'Year', y = 'Value', color = 'Model', fill = 'Model') +
    ggplot2::coord_cartesian(ylim = c(0,NA)) +
    theme_sablefish()

  return(ts_plot)
}

#' Get Fishery and Survey Selectivity Plots
#'
#' @param rep List of n_models of `SPoRC` report lists
#' @param model_names Vector of model names
#'
#' @returns Plots of terminal year fishery and survey selectivity by fleet, region, and sex across models
#' @export get_selex_plot
#'
#' @examples
#' \dontrun{
#' get_selex_plot(list(rep1, rep2), c("Model1", "Model2"))
#' }
get_selex_plot <- function(rep, model_names) {

  fishsel_plot_df <- data.frame()
  srvsel_plot_df <- data.frame()

  for(i in 1:length(rep)) {
    # Fishery Selectivity
    fishsel_plot_tmp_df <- reshape2::melt(rep[[i]]$fish_sel) %>%
      dplyr::rename(Region = Var1, Year = Var2, Age = Var3, Sex = Var4, Fleet = Var5) %>%
      dplyr::mutate(Region = paste("Region", Region),
                    Fleet = paste("Fleet", Fleet),
                    Sex = paste("Sex", Sex),
                    Model = model_names[i]
      )

    # Survey Selectivity
    srvsel_plot_tmp_df <- reshape2::melt(rep[[i]]$srv_sel) %>%
      dplyr::rename(Region = Var1, Year = Var2, Age = Var3, Sex = Var4, Fleet = Var5) %>%
      dplyr::mutate(Region = paste("Region", Region),
                    Fleet = paste("Fleet", Fleet),
                    Sex = paste("Sex", Sex),
                    Model = model_names[i]
      )

    fishsel_plot_df <- rbind(fishsel_plot_df, fishsel_plot_tmp_df)
    srvsel_plot_df <- rbind(srvsel_plot_df, srvsel_plot_tmp_df)

  }

  # fishery selectivity plot
  fish_sel_plot <- ggplot2::ggplot(fishsel_plot_df %>%
                                     dplyr::filter(Year == max(fishsel_plot_df$Year)),
                                   ggplot2::aes(x = Age, y = value, color = factor(Model))) +
    ggplot2::geom_line(lwd = 1.3) +
    ggplot2::facet_grid(Sex ~ Region + Fleet) +
    ggplot2::labs(x = 'Age', y = 'Terminal Fishery selectivity', color = 'Model') +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    theme_sablefish() +
    ggplot2::theme(legend.key.width = unit(2, "lines"))

  # survey selectivity plot
  srv_sel_plot <- ggplot2::ggplot(srvsel_plot_df %>%
                                     dplyr::filter(Year == max(srvsel_plot_df$Year)),
                                   ggplot2::aes(x = Age, y = value, color = factor(Model))) +
    ggplot2::geom_line(lwd = 1.3) +
    ggplot2::facet_grid(Sex ~ Region + Fleet) +
    ggplot2::labs(x = 'Age', y = 'Terminal Survey selectivity', color = 'Model') +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    theme_sablefish() +
    ggplot2::theme(legend.key.width = unit(2, "lines"))


  return(list(fish_sel_plot, srv_sel_plot))
}

#' Get Plots of Biological Quantities
#'
#' @param data List of n_models of `SPoRC` data lists
#' @param rep List of n_models of `SPoRC` report lists
#' @param model_names Vector of model names
#'
#' @returns A list of plots for terminal year movement, natural mortality, weight-at-age, and maturity at age across models
#' @export get_biological_plot
#'
#' @examples
#' \dontrun{
#' get_biological_plot(list(data1, data2), list(rep1, rep2), c("Model1", "Model2"))
#' }
get_biological_plot <- function(data,
                                rep,
                                model_names) {

  move_plot_df <- data.frame()
  natmort_plot_df <- data.frame()
  waa_plot_df <- data.frame()
  mataa_plot_df <- data.frame()

  for(i in 1:length(rep)) {

    # Movement
    move_plot_tmp_df <- reshape2::melt(rep[[i]]$Movement) %>%
      dplyr::rename(Region_From = Var1, Region_To = Var2, Year = Var3, Age = Var4, Sex = Var5) %>%
      dplyr::mutate(Region_From = paste("From Region", Region_From),
                    Region_To = paste("To Region", Region_To),
                    Sex = paste("Sex", Sex),
                    Model = model_names[i]
      )

    # Natural Mortality
    natmort_plot_tmp_df <- reshape2::melt(rep[[i]]$natmort) %>%
      dplyr::rename(Region = Var1, Year = Var2, Age = Var3, Sex = Var4) %>%
      dplyr::mutate(Region = paste("Region", Region),
                    Sex = paste("Sex", Sex),
                    Model = model_names[i]
      )

    # Weight-at-age
    waa_plot_tmp_df <- reshape2::melt(data[[i]]$WAA) %>%
      dplyr::rename(Region = Var1, Year = Var2, Age = Var3, Sex = Var4) %>%
      dplyr::mutate(Region = paste("Region", Region),
                    Sex = paste("Sex", Sex),
                    Model = model_names[i]
      )

    # Maturity at age
    mataa_plot_tmp_df <- reshape2::melt(data[[i]]$MatAA) %>%
      dplyr::rename(Region = Var1, Year = Var2, Age = Var3, Sex = Var4) %>%
      dplyr::mutate(Region = paste("Region", Region),
                    Sex = paste("Sex", Sex),
                    Model = model_names[i]
      )

    # bind all
    move_plot_df <- rbind(move_plot_df, move_plot_tmp_df)
    natmort_plot_df <- rbind(natmort_plot_df, natmort_plot_tmp_df)
    waa_plot_df <- rbind(waa_plot_df, waa_plot_tmp_df)
    mataa_plot_df <- rbind(mataa_plot_df, mataa_plot_tmp_df)
  }

  # Movement plot
  move_plot <- ggplot(move_plot_df %>% dplyr::filter(Year == max(move_plot_df$Year)),
                      ggplot2::aes(x = Age, y = value, color = factor(Model), lty = Sex)) +
    ggplot2::geom_line(lwd = 1) +
    ggplot2::facet_grid(Region_To~Region_From) +
    ggplot2::labs(x = 'Age', y = 'Movement Probabilities', color = 'Model', lty = 'Sex') +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    theme_sablefish() +
    ggplot2::theme(legend.key.width = unit(2, "lines"))

  # Natural mortality plot
  natmort_plot <- ggplot2::ggplot(natmort_plot_df %>%
                                    dplyr::filter(Year == max(natmort_plot_df$Year)),
                                  ggplot2::aes(x = Age, y = value, color = factor(Model))) +
    ggplot2::geom_line(lwd = 2) +
    ggplot2::facet_grid(Region~Sex) +
    ggplot2::labs(x = 'Age', y = 'Natural Mortality', color = 'Model') +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    theme_sablefish() +
    ggplot2::theme(legend.key.width = unit(2, "lines"))

  # Weight at age plot
  waa_plot <- ggplot2::ggplot(waa_plot_df %>%
                                dplyr::filter(Year == max(waa_plot_df$Year)),
                              ggplot2::aes(x = Age, y = value, color = factor(Model))) +
    ggplot2::geom_line(lwd = 2) +
    ggplot2::facet_grid(Region~Sex) +
    ggplot2::labs(x = 'Age', y = 'Weight at Age', color = 'Model') +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    theme_sablefish() +
    ggplot2::theme(legend.key.width = unit(2, "lines"))

  # Maturity plot
  mataa_plot <- ggplot2::ggplot(mataa_plot_df %>%
                                  dplyr::filter(Year == max(mataa_plot_df$Year)),
                                ggplot2::aes(x = Age, y = value, color = factor(Model))) +
    ggplot2::geom_line(lwd = 2) +
    ggplot2::facet_grid(Region~Sex) +
    ggplot2::labs(x = 'Age', y = 'Maturity at Age', color = 'Model') +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    theme_sablefish() +
    ggplot2::theme(legend.key.width = unit(2, "lines"))

  return(list(move_plot, natmort_plot, waa_plot, mataa_plot))

}

#' Get Data Fitted to Plot
#'
#' @param data List of n_models of `SPoRC` data lists
#' @param model_names Character vector of model names
#'
#' @returns A plot of data that were fitted to across models
#' @export get_data_fitted_plot
#'
#' @examples
#' \dontrun{
#' get_data_fitted_plot(list(data1, data2), c("Model1", "Model2"))
#' }
get_data_fitted_plot <- function(data,
                                 model_names
                                 ) {

  data_plot_all_df <- data.frame()
  for(i in 1:length(data)) {

    # Get tag release indicator
    if(data[[i]]$UseTagging == 1) {
      use_tag_indicator <- array(0, dim = c(max(data[[i]]$tag_release_indicator[,1]), max(data[[i]]$tag_release_indicator[,2])))
      use_tag_indicator[data[[i]]$tag_release_indicator[,1],data[[i]]$tag_release_indicator[,2]] <- 1
    }

    # Bind all data indicators together
    data_plot_df <- reshape2::melt(data[[i]]$UseSrvLenComps) %>% # Survey lengths
      dplyr::rename(Region = Var1, Year = Var2, Fleet = Var3) %>%
      dplyr::mutate(Type = paste('Survey Lengths', "Fleet", Fleet)) %>%

      dplyr::bind_rows(
        # survey ages
        reshape2::melt(data[[i]]$UseSrvAgeComps) %>%
          dplyr::rename(Region = Var1, Year = Var2, Fleet = Var3) %>%
          dplyr::mutate(Type = paste('Survey Ages', "Fleet", Fleet)),

        # fishery lengths
        reshape2::melt(data[[i]]$UseFishLenComps) %>%
          dplyr::rename(Region = Var1, Year = Var2, Fleet = Var3) %>%
          dplyr::mutate(Type = paste('Fishery Lengths', "Fleet", Fleet)),

        # fishery ages
        reshape2::melt(data[[i]]$UseFishAgeComps) %>%
          dplyr::rename(Region = Var1, Year = Var2, Fleet = Var3) %>%
          dplyr::mutate(Type = paste('Fishery Ages', "Fleet", Fleet)),

        # fishery catches
        reshape2::melt(data[[i]]$UseCatch) %>%
          dplyr::rename(Region = Var1, Year = Var2, Fleet = Var3) %>%
          dplyr::mutate(Type = paste('Fishery Catch', "Fleet", Fleet)),

        # fishery indices
        reshape2::melt(data[[i]]$UseFishIdx) %>%
          dplyr::rename(Region = Var1, Year = Var2, Fleet = Var3) %>%
          dplyr::mutate(Type = paste('Fishery Index', "Fleet", Fleet)),

        # survey indices
        reshape2::melt(data[[i]]$UseSrvIdx) %>%
          dplyr::rename(Region = Var1, Year = Var2, Fleet = Var3) %>%
          dplyr::mutate(Type = paste('Survey Index', "Fleet", Fleet))
      )

    # Add tagging if used
    if (data[[i]]$UseTagging == 1) {
      tag_df <- reshape2::melt(use_tag_indicator) %>%
        dplyr::rename(Region = Var1, Year = Var2) %>%
        dplyr::mutate(Type = 'Tagging', Fleet = NA)
      data_plot_df <- dplyr::bind_rows(data_plot_df, tag_df)
    }

    # remove data not fitted to
    data_plot_df <- data_plot_df %>%
      dplyr::filter(value != 0) %>%
      dplyr::mutate(Region = paste("Region", Region),
                    Model = model_names[i])

    data_plot_all_df <- rbind(data_plot_all_df, data_plot_df)
  }

  data_plot <- ggplot2::ggplot(data_plot_all_df,
                               ggplot2::aes(x = Year, y = Type, fill = Type)) +
    ggplot2::geom_point(size = 3, pch = 21, color = 'black', alpha = 0.8) +
    ggplot2::facet_grid(Model~Region) +
    theme_sablefish() +
    ggplot2::theme(legend.position = 'none') +
    ggplot2::labs(x = 'Year', y = '')

  return(data_plot)
}

#' Get plot of negative log likelihood values
#'
#' @param rep List of n_models of `SPoRC` report lists
#' @param model_names Vector of model names
#'
#' @returns Plot and tables of negative log likelihood values across models
#' @export get_nLL_plot
#'
#' @examples
#' \dontrun{
#' get_nLL_plot(list(rep1, rep2), c("Model1", "Model2"))
#' }
get_nLL_plot <- function(rep,
                         model_names
                         ) {


  nLL_all_df <- data.frame() # empty dataframe
  for(i in 1:length(rep)) {

    # Negative log likelihoods
    nLL_df <- data.frame(

      # nLL values
      value = c(rep[[i]]$jnLL,
                rep[[i]]$h_nLL,
                rep[[i]]$M_nLL,
                sum(rep[[i]]$Rec_nLL),
                rep[[i]]$sel_nLL,
                sum(rep[[i]]$Tag_nLL),
                sum(rep[[i]]$Catch_nLL),
                sum(rep[[i]]$Fmort_nLL),
                rep[[i]]$srv_q_nLL,
                rep[[i]]$fish_q_nLL,
                sum(rep[[i]]$SrvIdx_nLL),
                rep[[i]]$TagRep_nLL,
                sum(rep[[i]]$FishIdx_nLL),
                sum(rep[[i]]$Init_Rec_nLL),
                rep[[i]]$Movement_nLL,
                sum(rep[[i]]$SrvAgeComps_nLL),
                sum(rep[[i]]$FishAgeComps_nLL),
                sum(rep[[i]]$SrvLenComps_nLL),
                sum(rep[[i]]$FishLenComps_nLL)),

      # nLL names
      name = c("jnLL", "Steepness Prior", "M Prior", "Recruitment Penalty",
               "Selectivity Penalty", "Tag nLL", "Catch nLL", "Fishing Mortality Penalty",
               "Survey Q Prior", "Fishery Q Prior", "Survey Index nLL", "Tag Reporting Prior",
               "Fishery Index nLL", "Initial Age Penalty", "Movement Prior",
               "Survey Age nLL", "Fishery Age nLL", "Survey Length nLL", "Fishery Length nLL"),

      # nLL types
      type = c('jnLL', 'Prior', "Prior", "Penalty", "Penalty", "Tagging",
               "Catch", "Penalty", "Prior", "Prior", "Index",
               "Prior", "Index", "Penalty", "Prior", "Age",
               "Age", "Length", "Length"),

      Model = model_names[i]
    )

    nLL_all_df <- rbind(nLL_all_df, nLL_df)
  }

  nLL_plot <- ggplot2::ggplot(nLL_all_df %>% dplyr::filter(value != 0),
                              ggplot2::aes(x = name, y = value, fill = type)) +
    ggplot2::geom_col() +
    ggplot2::facet_wrap(~Model) +
    theme_sablefish() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.5)) +
    ggplot2::labs(x = 'Likelihood Component', y = 'Likelihood', fill = 'Type')


  # get nLL table
  nLL_table <- nLL_all_df %>%
    dplyr::select(-type) %>%
    dplyr::filter(value != 0) %>%  # remove zeros
    tidyr::pivot_wider(names_from = name, values_from = value)
  table_plot <- grid::grid.grabExpr(gridExtra::grid.table(nLL_table))
  table_plot1 <- cowplot::ggdraw() + cowplot::draw_grob(table_plot)

  return(list(nLL_plot, table_plot1))
}

#' Get Index Fits Plot
#'
#' @param data List of n_models of `SPoRC` data lists
#' @param rep List of n_models of `SPoRC` report lists
#' @param model_names Vector of model names
#'
#' @returns A plot of fitted values to various indices across models
#' @export get_idx_fits_plot
#'
#' @examples
#' \dontrun{
#' get_idx_fits_plot(list(data1, data2), list(rep1, rep2), c("Model1", "Model2"))
#' }
get_idx_fits_plot <- function(data,
                              rep,
                              model_names
                              ) {

  idx_fits_all <- data.frame()
  # get index fits data
  for(i in 1:length(rep)) {
    idx_fits <- get_idx_fits(data = data[[i]], rep = rep[[i]], year_labs = data[[i]]$years) %>%
      dplyr::mutate(Model = model_names[i])
    idx_fits_all <- rbind(idx_fits_all, idx_fits) # bind
  }

  # Plot index fits
  idx_fit_plot <- ggplot2::ggplot() +
    ggplot2::geom_line(idx_fits_all %>% dplyr::filter(obs != 0),
                       mapping = ggplot2::aes(x = Year, y = value, color = factor(Model)), lwd = 1.3) +
    ggplot2::geom_pointrange(idx_fits_all %>% dplyr::filter(obs != 0),
                             mapping = ggplot2::aes(x = Year, y = obs, ymin = lci, ymax = uci), color = 'black') +
    ggplot2::labs(x = "Year", y = 'Index', color = 'Model') +
    theme_sablefish() +
    ggplot2::coord_cartesian(ylim = c(0,NA)) +
    ggplot2::facet_grid(Category~Region, scales = 'free_y')

  return(idx_fit_plot)
}

#' Get Retrospective Plot
#'
#' @param retro_output Dataframe generated from do_retrospective
#' @param Rec_Age Age in which recruitment occurs
#'
#' @returns A retrospective plot of recruitment and SSB in relative and absolute scales, as well as a retrospective plot of recruitment by cohort (squid plot)
#' @export get_retrospective_plot
#'
#' @examples
#' \dontrun{
#' # do retrospective
#' retro <- do_retrospective(n_retro = 7, # number of retro peels to run
#' data = data, # rtmb data
#' parameters = parameters, # rtmb parameters
#' mapping = mapping, # rtmb mapping
#' random = NULL, # if random effects are used
#' do_par = TRUE, # whether or not to parralleize
#' n_cores = 7, # if parallel, number of cores to use
#' do_francis = F, # if we want tod o Francis
#' n_francis_iter = NULL # Number of francis iterations to do
#' )
#' get_retrospective_plot(retro, Rec_Age = 2)
#' }
get_retrospective_plot <- function(retro_output, Rec_Age) {

  if(sum(names(retro_output) %in% c("Region", "Year", "value", "Type", "peel")) != 5) stop("Retro output does not contain Region, Year, value, Type, peel.")

  # Get relative differences
  ret_df <- get_retrospective_relative_difference(retro_output)

  # Get mohns rho (mean relative difference for a given terminal year peel to terminal year estimates)
  mohns_rho <- ret_df %>%
    dplyr::filter(peel == (max(Year) - Year)) %>%
    dplyr::group_by(Type, Region) %>%
    dplyr::summarize(rho = mean(rd))

  # get retrospective plot
  retro_plot <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, lty = 2, lwd = 1.3) +
    ggplot2::geom_line(ret_df,
                       mapping = ggplot2::aes(x = Year, y = rd, group = as.numeric(peel), color = as.numeric(peel)), lwd = 1.3) +
    ggplot2::geom_point(ret_df %>% dplyr::filter(peel == max(Year) - Year),
                        mapping = ggplot2::aes(x = Year, y = rd, group = as.numeric(peel), fill = as.numeric(peel)),
                        pch = 21, size = 6) +
    ggplot2::geom_text(mohns_rho, mapping = aes(x = -Inf, y = Inf, label = paste("Mohns Rho:", round(rho, 4))),
                       hjust = -0.3, vjust = 3, size = 5) +
    ggplot2::guides(color = ggplot2::guide_colourbar(barwidth = 15, barheight = 1.3)) +
    ggplot2::labs(x = 'Year', y = 'Relative Difference from Terminal Year', color = 'Retrospective Year', fill = 'Retrospective Year') +
    ggplot2::scale_color_viridis_c() +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::facet_grid(Region~Type, scales = 'free') +
    theme_sablefish() +
    ggplot2::theme(legend.position = 'top')

  # get absolute retro plot
  abs_retro_plot <- ggplot2::ggplot() +
    ggplot2::geom_line(retro_output %>% filter(peel != 0), mapping = ggplot2::aes(x = Year, y = value, group = peel, color = peel), lwd = 1) +
    ggplot2::geom_line(retro_output %>% filter(peel == 0), mapping = ggplot2::aes(x = Year, y = value), lty = 2, lwd = 1) +
    ggplot2::scale_color_viridis_c() +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    ggplot2::facet_wrap(~Type, scales = 'free') +
    theme_sablefish() +
    ggplot2::labs(x = 'Year', y = 'Value', color = 'Peel')

  # get squid plot
  squid_plot <- retro %>%
    dplyr::mutate(Year = Year, terminal = max(retro$Year) - peel, cohort = Year - Rec_Age, years_est = terminal-Year) %>%
    dplyr::filter(Type == 'Recruitment', cohort %in% seq(max(retro$Year) - 10, max(retro$Year), 1), terminal != Year) %>%
    ggplot2::ggplot(ggplot2::aes(x = years_est - 1, y = value, group = Year, color = factor(cohort))) +
    ggplot2::geom_line(lwd = 1.3) +
    ggplot2::geom_point(size = 4) +
    ggplot2::theme_bw(base_size = 15) +
    ggplot2::labs(x = 'Years since cohort was last estimated', y = 'Recruitment (millions)', color = 'Cohort')

  return(list(retro_plot, abs_retro_plot, squid_plot))
}

#' Plotting function for all basic quantities
#'
#' @param data List of n_models of `SPoRC` data lists
#' @param rep List of n_models of `SPoRC` report lists
#' @param sd_rep List of n_models of sd report lists from `SPoRC`
#' @param out_path Path to the output directory. Users only need to specify the path.
#' @param model_names Character vector of model names
#'
#' @returns A series of plots compared across models outputted as a pdf in the specified directory
#' @export plot_all_basic
#'
#' @examples
#' \dontrun{
#' plot_all_basic(
#'   data = list(data1, data2),
#'   rep = list(rep1, rep2),
#'   sd_rep = list(sd_rep1, sd_rep2),
#'   model_names = c("Model1", "Model2"),
#'   out_path = here::here()
#' )
#' }
plot_all_basic <- function(data,
                           rep,
                           sd_rep,
                           model_names,
                           out_path) {

  pdf(here::here(out_path, "plot_results.pdf"), width = 25, height = 13)
  print(get_biological_plot(data = data, rep = rep, model_names = model_names))
  print(get_data_fitted_plot(data = data, model_names = model_names))
  print(get_ts_plot(rep = rep, sd_rep = sd_rep, model_names = model_names))
  print(get_selex_plot(rep = rep, model_names = model_names))
  print(get_nLL_plot(rep = rep, model_names = model_names))
  dev.off()

}

