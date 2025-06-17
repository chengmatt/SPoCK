#' Title Get Time Series Plots
#'
#' @param rep RTMB report list
#' @param sd_rep RTMB sdreport list
#'
#' @returns Plots of spawning biomass, total biomass, recruitment, and fishing mortality time-series
#' @export get_ts_plot
#'
#' @examples
#' \dontrun{
#'   get_ts_plot(rep, sd_rep)
#' }
get_ts_plot <- function(rep,
                         sd_rep
                         ) {

  # Spawning Stock Biomass
  ssb_plot_df <- reshape2::melt(rep$SSB) %>%
    dplyr::rename(Region = Var1, Year = Var2) %>%
    dplyr::bind_cols(se = sd_rep$sd[names(sd_rep$value) == "log(SSB)"]) %>%
    dplyr::mutate(lwr = exp(log(value) - 1.96 * se),
                  upr = exp(log(value) + 1.96 * se),
                  Region = paste("Region", Region),
                  Type = 'SSB')

  # Total Biomass
  totbiom_plot_df <- reshape2::melt(rep$Total_Biom) %>%
    dplyr::rename(Region = Var1, Year = Var2) %>%
    dplyr::bind_cols(se = sd_rep$sd[names(sd_rep$value) == "log(Total_Biom)"]) %>%
    dplyr::mutate(lwr = exp(log(value) - 1.96 * se),
                  upr = exp(log(value) + 1.96 * se),
                  Region = paste("Region", Region),
                  Type = 'Total Biom')

  # Recruitment
  rec_plot_df <- reshape2::melt(rep$Rec) %>%
    dplyr::rename(Region = Var1, Year = Var2) %>%
    dplyr::bind_cols(se = sd_rep$sd[names(sd_rep$value) == "log(Rec)"]) %>%
    dplyr::mutate(lwr = exp(log(value) - 1.96 * se),
                  upr = exp(log(value) + 1.96 * se),
                  Region = paste("Region", Region),
                  Type = 'Recruitment')

  # Fishing Mortality
  f_plot_df <- reshape2::melt(rep$Fmort) %>%
    dplyr::rename(Region = Var1, Year = Var2, Type = Var3) %>%
    dplyr::mutate(Region = paste("Region", Region),
                  Type = paste("Fleet", Type, "F"),
                  lwr = NA,
                  upr = NA,
                  se = NA)

  # bind together
  biom_rec_df <- rbind(ssb_plot_df, totbiom_plot_df, rec_plot_df, f_plot_df)

  # Plot time series
  ts_plot <- ggplot2::ggplot(biom_rec_df,
                             ggplot2::aes(x = Year, y = value, ymin = lwr, ymax = upr)) +
    ggplot2::geom_line(lwd = 0.9) +
    ggplot2::geom_ribbon(alpha = 0.3) +
    ggplot2::facet_grid(Type~Region, scales = 'free') +
    ggplot2::labs(x = 'Year', y = 'Value') +
    ggplot2::coord_cartesian(ylim = c(0,NA)) +
    theme_sablefish()

  return(ts_plot)
}

#' Title Get Fishery and Survey Selectivity Plots
#'
#' @param rep RTMB report list
#'
#' @returns Plots of fishery and survey selectivity by fleet, region, and sex
#' @export get_selex_plot
#'
#' @examples
#' \dontrun{
#' get_selex_plot(rep)
#' }
get_selex_plot <- function(rep) {

  # Fishery Selectivity
  fishsel_plot_df <- reshape2::melt(rep$fish_sel) %>%
    dplyr::rename(Region = Var1, Year = Var2, Age = Var3, Sex = Var4, Fleet = Var5) %>%
    dplyr::mutate(Region = paste("Region", Region),
                  Fleet = paste("Fleet", Fleet),
                  Sex = paste("Sex", Sex)
    )

  # Survey Selectivity
  srvsel_plot_df <- reshape2::melt(rep$srv_sel) %>%
    dplyr::rename(Region = Var1, Year = Var2, Age = Var3, Sex = Var4, Fleet = Var5) %>%
    dplyr::mutate(Region = paste("Region", Region),
                  Fleet = paste("Fleet", Fleet),
                  Sex = paste("Sex", Sex)
    )

  # fishery selectivity plot
  fish_sel_plot <- ggplot2::ggplot(fishsel_plot_df %>%
                                     dplyr::filter(Year %in% c(seq(min(fishsel_plot_df$Year), max(fishsel_plot_df$Year), 5))),
                                   ggplot2::aes(x = Age, y = value, color = Year, group = Year)) +
    ggplot2::geom_line(lwd = 1.3) +
    ggplot2::facet_grid(Sex ~ Region + Fleet) +
    ggplot2::scale_color_viridis_c() +
    ggplot2::labs(x = 'Age', y = 'Fishery selectivity') +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    theme_sablefish() +
    ggplot2::theme(legend.key.width = unit(2, "lines"))

  # survey selectivity plot
  srv_sel_plot <- ggplot2::ggplot(srvsel_plot_df %>%
                                    dplyr::filter(Year %in% c(seq(min(srvsel_plot_df$Year), max(srvsel_plot_df$Year), 5))),
                                  ggplot2::aes(x = Age, y = value, color = Year, group = Year)) +
    ggplot2::geom_line(lwd = 1.3) +
    ggplot2::facet_grid(Sex ~ Region + Fleet) +
    ggplot2::scale_color_viridis_c() +
    ggplot2::labs(x = 'Age', y = 'Survey selectivity') +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    theme_sablefish() +
    ggplot2::theme(legend.key.width = unit(2, "lines"))


  return(list(fish_sel_plot, srv_sel_plot))
}

#' Title Get Plots of Biological Quantities
#'
#' @param data data list used for SPoCK
#' @param rep report list used for SPoCK
#'
#' @returns A list of plots for movement, natural mortality, weight-at-age, and maturity at age
#' @export get_biological_plot
#'
#' @examples
#' \dontrun{
#' get_biological_plot(data, rep)
#' }
get_biological_plot <- function(data,
                                 rep) {

  # Movement
  move_plot_df <- reshape2::melt(rep$Movement) %>%
    dplyr::rename(Region_From = Var1, Region_To = Var2, Year = Var3, Age = Var4, Sex = Var5) %>%
    dplyr::mutate(Region_From = paste("From Region", Region_From),
                  Region_To = paste("To Region", Region_To),
                  Sex = paste("Sex", Sex)
    )

  # Natural Mortality
  natmort_plot_df <- reshape2::melt(rep$natmort) %>%
    dplyr::rename(Region = Var1, Year = Var2, Age = Var3, Sex = Var4) %>%
    dplyr::mutate(Region = paste("Region", Region),
                  Sex = paste("Sex", Sex)
    )

  # Weight-at-age
  waa_plot_df <- reshape2::melt(data$WAA) %>%
    dplyr::rename(Region = Var1, Year = Var2, Age = Var3, Sex = Var4) %>%
    dplyr::mutate(Region = paste("Region", Region),
                  Sex = paste("Sex", Sex)
    )

  # Maturity at age
  mataa_plot_df <- reshape2::melt(data$MatAA) %>%
    dplyr::rename(Region = Var1, Year = Var2, Age = Var3, Sex = Var4) %>%
    dplyr::mutate(Region = paste("Region", Region),
                  Sex = paste("Sex", Sex)
    )

  # Movement plot
  move_plot <- ggplot(move_plot_df,
                      ggplot2::aes(x = Year, y = Age, fill = value, z = value)) +
    ggplot2::geom_tile(lwd = 1) +
    ggplot2::facet_grid(Region_To~Region_From + Sex) +
    ggplot2::labs(x = 'Year', y = 'Age', fill = 'Movement Probabilities') +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    theme_sablefish() +
    scale_fill_viridis_c(option = 'magma') +
    ggplot2::theme(legend.key.width = unit(2, "lines"))

  # Natural mortality plot
  natmort_plot <- ggplot2::ggplot(natmort_plot_df %>%
                                    dplyr::filter(Year %in% c(seq(min(natmort_plot_df$Year), max(natmort_plot_df$Year), 5))),
                                  ggplot2::aes(x = Age, y = value, color = Year)) +
    ggplot2::geom_line(lwd = 2) +
    ggplot2::facet_grid(Region~Sex) +
    ggplot2::scale_color_viridis_c() +
    ggplot2::labs(x = 'Age', y = 'Natural Mortality', color = 'Year') +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    theme_sablefish() +
    ggplot2::theme(legend.key.width = unit(2, "lines"))

  # Weight at age plot
  waa_plot <- ggplot2::ggplot(waa_plot_df %>%
                                dplyr::filter(Year %in% c(seq(min(waa_plot_df$Year), max(waa_plot_df$Year), 5))),
                              ggplot2::aes(x = Age, y = value, color = Year, group = Year)) +
    ggplot2::geom_line(lwd = 2) +
    ggplot2::facet_grid(Region~Sex) +
    ggplot2::scale_color_viridis_c() +
    ggplot2::labs(x = 'Age', y = 'Weight at Age', color = 'Year') +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    theme_sablefish() +
    ggplot2::theme(legend.key.width = unit(2, "lines"))

  # Maturity plot
  mataa_plot <- ggplot2::ggplot(mataa_plot_df %>%
                                  dplyr::filter(Year %in% c(seq(min(mataa_plot_df$Year), max(mataa_plot_df$Year), 5))),
                                ggplot2::aes(x = Age, y = value, color = Year, group = Year)) +
    ggplot2::geom_line(lwd = 2) +
    ggplot2::facet_grid(Region~Sex) +
    ggplot2::scale_color_viridis_c() +
    ggplot2::labs(x = 'Age', y = 'Maturity at Age', color = 'Year') +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    theme_sablefish() +
    ggplot2::theme(legend.key.width = unit(2, "lines"))

  return(list(move_plot, natmort_plot, waa_plot, mataa_plot))

}

#' Title Get Data Fitted to Plot
#'
#' @param data data list fed into SPoCK
#'
#' @returns A plot of data that were fitted to
#' @export get_data_fitted_plot
#'
#' @examples
#' \dontrun{
#' get_data_fitted_plot(data)
#' }
get_data_fitted_plot <- function(data) {

  # Get tag release indicator
  if(data$UseTagging == 1) {
    use_tag_indicator <- array(0, dim = c(max(data$tag_release_indicator[,1]), max(data$tag_release_indicator[,2])))
    use_tag_indicator[data$tag_release_indicator[,1],data$tag_release_indicator[,2]] <- 1
  }

  # Bind all data indicators together
  data_plot_df <- reshape2::melt(data$UseSrvLenComps) %>% # Survey lengths
    dplyr::rename(Region = Var1, Year = Var2, Fleet = Var3) %>%
    dplyr::mutate(Type = paste('Survey Lengths', "Fleet", Fleet)) %>%

    dplyr::bind_rows(
      # survey ages
      reshape2::melt(data$UseSrvAgeComps) %>%
        dplyr::rename(Region = Var1, Year = Var2, Fleet = Var3) %>%
        dplyr::mutate(Type = paste('Survey Ages', "Fleet", Fleet)),

      # fishery lengths
      reshape2::melt(data$UseFishLenComps) %>%
        dplyr::rename(Region = Var1, Year = Var2, Fleet = Var3) %>%
        dplyr::mutate(Type = paste('Fishery Lengths', "Fleet", Fleet)),

      # fishery ages
      reshape2::melt(data$UseFishAgeComps) %>%
        dplyr::rename(Region = Var1, Year = Var2, Fleet = Var3) %>%
        dplyr::mutate(Type = paste('Fishery Ages', "Fleet", Fleet)),

      # fishery catches
      reshape2::melt(data$UseCatch) %>%
        dplyr::rename(Region = Var1, Year = Var2, Fleet = Var3) %>%
        dplyr::mutate(Type = paste('Fishery Catch', "Fleet", Fleet)),

      # fishery indices
      reshape2::melt(data$UseFishIdx) %>%
        dplyr::rename(Region = Var1, Year = Var2, Fleet = Var3) %>%
        dplyr::mutate(Type = paste('Fishery Index', "Fleet", Fleet)),

      # survey indices
      reshape2::melt(data$UseSrvIdx) %>%
        dplyr::rename(Region = Var1, Year = Var2, Fleet = Var3) %>%
        dplyr::mutate(Type = paste('Survey Index', "Fleet", Fleet))
    )

  # Add tagging if used
  if (data$UseTagging == 1) {
    tag_df <- reshape2::melt(use_tag_indicator) %>%
      dplyr::rename(Region = Var1, Year = Var2) %>%
      dplyr::mutate(Type = 'Tagging', Fleet = NA)
    data_plot_df <- dplyr::bind_rows(data_plot_df, tag_df)
  }

  # remove data not fitted to
  data_plot_df <- data_plot_df %>%
    dplyr::filter(value != 0) %>%
    dplyr::mutate(Region = paste("Region", Region))

  data_plot <- ggplot2::ggplot(data_plot_df,
                               ggplot2::aes(x = Year, y = Type, fill = Type)) +
    ggplot2::geom_point(size = 3, pch = 21, color = 'black', alpha = 0.8) +
    ggplot2::facet_wrap(~Region) +
    theme_sablefish() +
    ggplot2::theme(legend.position = 'none') +
    ggplot2::labs(x = 'Year', y = '')

  return(data_plot)
}

#' Title Get plot of negative log likelihood values
#'
#' @param rep Report list from SPoCK
#' @param sd_rep sdreport list from SPoCK
#'
#' @returns Plot of negative log likelihood values
#' @export get_nLL_plot
#'
#' @examples
#' \dontrun{
#' get_nLL_plot(rep, sd_rep)
#' }
get_nLL_plot <- function(rep, sd_rep) {

  # Negative log likelihoods
  nLL_df <- data.frame(

    # nLL values
    value = c(rep$jnLL,
              rep$h_nLL,
              rep$M_nLL,
              sum(rep$Rec_nLL),
              rep$sel_nLL,
              sum(rep$Tag_nLL),
              sum(rep$Catch_nLL),
              sum(rep$Fmort_nLL),
              rep$srv_q_nLL,
              rep$fish_q_nLL,
              sum(rep$SrvIdx_nLL),
              rep$TagRep_nLL,
              sum(rep$FishIdx_nLL),
              sum(rep$Init_Rec_nLL),
              rep$Movement_nLL,
              sum(rep$SrvAgeComps_nLL),
              sum(rep$FishAgeComps_nLL),
              sum(rep$SrvLenComps_nLL),
              sum(rep$FishLenComps_nLL)),

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
             "Age", "Length", "Length")
  )

  nLL_plot <- ggplot2::ggplot(nLL_df %>% dplyr::filter(value != 0),
                              ggplot2::aes(x = name, y = value, fill = type)) +
    ggplot2::geom_col() +

    # parameters estimated
    ggplot2::annotate('text', x = -Inf, y = -Inf,
                      label = paste("Parameters Est:", length(sd_rep$par.fixed)), size = 8,
                      hjust = -0.17, vjust = -11) +
    # gradient
    ggplot2::annotate('text', x = -Inf, y = -Inf,
                      label = paste("Gradient:", round(max(sd_rep$gradient.fixed), 10)), size = 8,
                      hjust = -0.33, vjust = -13) +

    # gradient
    ggplot2::annotate('text', x = -Inf, y = -Inf,
                      label = paste("Hessian:", sd_rep$pdHess), size = 8,
                      hjust = -0.25, vjust = -15) +

    theme_sablefish() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.5)) +
    ggplot2::labs(x = 'Likelihood Component', y = 'Likelihood', fill = 'Type')

  return(nLL_plot)
}

#' Title Get Index Fits Plot
#'
#' @param data data list from SPoCK
#' @param rep report list from SPoCK
#'
#' @returns A plot of fitted values to various indices
#' @export get_idx_fits_plot
#'
#' @examples
#' \dontrun{
#' get_idx_fits_plot(data,rep)
#' }
get_idx_fits_plot <- function(data,
                               rep) {

  # get index fits data
  idx_fits <- get_idx_fits(data = data, rep = rep, year_labs = data$years)

  # Plot index fits
  idx_fit_plot <- ggplot2::ggplot() +
    ggplot2::geom_line(idx_fits %>% dplyr::filter(obs != 0),
                       mapping = ggplot2::aes(x = Year, y = value), lwd = 1.3, col = 'black') +
    ggplot2::geom_pointrange(idx_fits %>% dplyr::filter(obs != 0),
                             mapping = ggplot2::aes(x = Year, y = obs, ymin = lci, ymax = uci), color = 'blue') +
    ggplot2::labs(x = "Year", y = 'Index') +
    theme_sablefish() +
    ggplot2::coord_cartesian(ylim = c(0,NA)) +
    ggplot2::facet_grid(Category~Region, scales = 'free_y')

  return(idx_fit_plot)
}

#' Title Get Retrospective Plot
#'
#' @param retro_output Dataframe generated from do_retrospective
#'
#' @returns A retrospective plot of recruitment and SSB
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
#' get_retrospective_plot(retro)
#' }
get_retrospective_plot <- function(retro_output) {

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

  return(retro_plot)
}

