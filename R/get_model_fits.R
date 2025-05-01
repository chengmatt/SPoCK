#' Gets index fits results
#'
#' @param data Data list fed into RTMB
#' @param year_labs Year labels to use (vector)
#' @param rep Report list output from RTMB
#'
#' @returns Fits to indices as a dataframe
#' @export get_idx_fits
#' @import dplyr
#' @importFrom tidyr drop_na
#' @examples
#' \dontrun{
#' idx_fits <- get_idx_fits(data = data, rep = rep, year_labs = seq(1960, 2024, 1))
#'
#' idx_fits <- idx_fits %>%
#'   mutate(
#'     Idx = case_when(
#'       Type == "Fishery" & Year < 1995 ~ "Japanese Fishery CPUE Index",
#'       Type == "Fishery" & Year >= 1995 ~ "Domestic Fishery CPUE Index",
#'       Type == 'Survey' & Fleet == 1 ~ "Domestic LL Survey Relative Population Numbers",
#'       Type == 'Survey' & Fleet == 2 ~ "GOA Trawl Survey Biomass (kt)",
#'       Type == 'Survey' & Fleet == 3 ~ 'Japanese LL Survey Relative Population Numbers'
#'     )
#'   )
#' ggplot() +
#'   geom_line(idx_fits, mapping = aes(x = Year, y = value), lwd = 1.3, col = 'red') +
#'   geom_pointrange(idx_fits, mapping = aes(x = Year, y = obs, ymin = lci, ymax = uci), color = 'blue', pch = 1) +
#'   labs(x = "Year", y = 'Index') +
#'   theme_bw(base_size = 20) +
#'   facet_wrap(~Idx, scales = 'free', ncol = 2)
#'
#' }
get_idx_fits <- function(data,
                         rep,
                         year_labs
                         ) {

  colnames(data$ObsSrvIdx) <- year_labs
  colnames(data$ObsSrvIdx_SE) <- year_labs
  colnames(rep$PredSrvIdx) <- year_labs
  colnames(data$ObsFishIdx) <- year_labs
  colnames(data$ObsFishIdx_SE) <- year_labs
  colnames(rep$PredFishIdx) <- year_labs

  # Observed survey index
  obs_srv <- reshape2::melt(data$ObsSrvIdx) %>% dplyr::rename(obs = value) %>%
    dplyr::left_join(reshape2::melt(data$ObsSrvIdx_SE) %>%  dplyr::rename(se = value), by = c("Var1", "Var2", "Var3")) %>%
    dplyr::mutate(lci = exp(log(obs) - (1.96 * se)), uci = exp(log(obs) + (1.96 * se)), Type = 'Survey') %>%
    tidyr::drop_na() %>%
    dplyr::rename(Region = Var1, Year = Var2, Fleet = Var3)

  # Predicted survey index
  pred_srv <- reshape2::melt(rep$PredSrvIdx) %>%
    dplyr::rename(Region = Var1, Year = Var2, Fleet = Var3) %>%
    dplyr::mutate(Type = 'Survey') %>%
    dplyr::filter(Year %in% unique(obs_srv$Year), value != 0)

  # combine survey results
  all_srv <- obs_srv %>%
    dplyr::left_join(pred_srv, by = c("Region", "Year", "Fleet", 'Type')) %>%
    dplyr::mutate(resid = log(obs) - log(value))

  # Observed fishery index
  obs_fish <- reshape2::melt(data$ObsFishIdx) %>%
    dplyr::rename(obs = value) %>%
    dplyr::left_join(reshape2::melt(data$ObsFishIdx_SE) %>%
                dplyr::rename(se = value), by = c("Var1", "Var2", "Var3")) %>%
    dplyr::mutate(lci = exp(log(obs) - (1.96 * se)), uci = exp(log(obs) + (1.96 * se)), Type = 'Fishery') %>%
    tidyr::drop_na() %>%
    dplyr::rename(Region = Var1, Year = Var2, Fleet = Var3)

  # Predicted fishery index
  pred_fish <- reshape2::melt(rep$PredFishIdx) %>%
    dplyr::rename(Region = Var1, Year = Var2, Fleet = Var3) %>%
    dplyr::mutate(Type = 'Fishery') %>%
    dplyr::filter(Year %in% unique(obs_fish$Year), value != 0)

  # combine survey results
  all_fish <- obs_fish %>% dplyr::left_join(pred_fish, by = c("Region", "Year", "Fleet", 'Type')) %>%
    dplyr::mutate(resid = log(obs) - log(value))

  all_idx <- rbind(all_fish, all_srv)
  return(all_idx)
}

#' Restructure composition values, used within a variety of functions to either do Francis reweighting or get observed and expected composition values
#'
#' @param Exp Expected values (catch at age or survey index at age) indexed for a given year and fleet (structured as a matrix by age and sex)
#' @param Obs Observed values (catch at age or survey index at age) indexed for a given year and fleet (structured as a matrix by age and sex)
#' @param Comp_Type Composition Parameterization Type (== 0, aggregated comps by sex, == 1, split comps by sex (no implicit sex ratio information), == 2, joint comps across sexes (implicit sex ratio information), == 3 joint comps across sexes and regions (implicit sex ratio and region information))
#' @param age_or_len Age or length comps (== 0, Age, == 1, Length)
#' @param AgeingError Ageing Error matrix
#' @param comp_agg_type Composition aggregation type to mimic sablefish ADMB assessment
#'
#' @return Returns a list of array observed and expected values for a given year and fleet
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' comps <- Restrc_Comps(Exp, Obs, Comp_Type, age_or_len, AgeingError)
#' comps$Exp; comps$Obs
#' }
Restrc_Comps <- function(Exp,
                         Obs,
                         Comp_Type,
                         age_or_len,
                         AgeingError,
                         comp_agg_type
                         ) {

  const <- 0

  # Add constant to observed and expected (gets normalized later on)
  Obs <- Obs + const
  Exp <- Exp + const

  # Dimensions
  n_regions <- dim(Exp)[1]
  n_bins <- dim(Exp)[3]
  n_sexes <- dim(Exp)[4]

  # Storage
  Exp_mat = array(NA, c(n_regions, n_bins, n_sexes))
  Obs_mat = array(NA, c(n_regions, n_bins, n_sexes))

  # Aggregated comps by sex and region
  if(Comp_Type == 0) {
    if(comp_agg_type == 0) { # aggregated age comps are normalized, aggregated, ageing error, and then normalized again
      # Expected Values
      tmp_Exp = Exp / array(data = rep(colSums(matrix(Exp, nrow = n_bins)), each = n_bins), dim = dim(Exp)) # normalize by sex and region
      tmp_Exp = matrix(rowSums(matrix(tmp_Exp, nrow = n_bins)) / (n_sexes * n_regions), nrow = 1) # take average proportions and transpose
    }

    if(comp_agg_type == 1) tmp_Exp = matrix(rowSums(matrix(Exp, nrow = n_bins)) / (n_sexes * n_regions), nrow = 1) # age comps are aggregated, ageing error, and the normalized

    if(age_or_len == 0) {
      tmp_Exp = tmp_Exp %*% AgeingError # apply ageing error
      tmp_Exp = as.vector((tmp_Exp) / sum(tmp_Exp)) # renormalize
    }
    if(age_or_len == 1) tmp_Exp = as.vector((tmp_Exp) / sum(tmp_Exp)) # renormalize (lengths)

    # Normalize observed
    tmp_Obs =  (Obs[1,1,,1,1]) / sum( Obs[1,1,,1,1])

    # Input into storage matrix
    Exp_mat[1,,1] = tmp_Exp
    Obs_mat[1,,1] = tmp_Obs

  } # end if aggregated comps across sex

  # 'Split' comps by sex (no implicit sex ratio information)
  if(Comp_Type == 1) {
    for(s in 1:n_sexes) {
      for(r in 1:n_regions) {

        # Expected Values
        if(age_or_len == 0) tmp_Exp = ((Exp[r,1,,s,1]) / sum(Exp[r,1,,s,1])) %*% AgeingError # Normalize temporary variable (ages)
        if(age_or_len == 1 && comp_agg_type == 0) tmp_Exp = (Exp[r,1,,s,1]) / sum(Exp[r,1,,s,1]) # Length comps are not normalized prior to age length transition
        if(age_or_len == 1 && comp_agg_type == 1) tmp_Exp = (Exp[r,1,,s,1]) # Length comps are normalized prior to age length transition

        tmp_Obs = (Obs[r,1,,s,1]) / sum(Obs[r,1,,s,1]) # Normalize observed temporary variable

        # Input into storage matrix
        Exp_mat[r,,s] = tmp_Exp
        Obs_mat[r,,s] = tmp_Obs
      } # end r loop
    } # end s loop
  } # end if 'Split' comps by sex

  if(Comp_Type == 2) {
    for(r in 1:n_regions) {
      # Expected values
      if(age_or_len == 0) { # if ages
        tmp_Exp = t(as.vector((Exp[r,1,,,1])/ sum(Exp[r,1,,,1]))) %*% kronecker(diag(n_sexes), AgeingError) # apply ageing error
        tmp_Exp = as.vector((tmp_Exp) / sum(tmp_Exp)) # renormalize to make sure sum to 1
      } # if ages

      if(age_or_len == 1) tmp_Exp = as.vector((Exp[r,1,,,1]) / sum((Exp[r,1,,,1]))) # Normalize temporary variable (lengths)

      tmp_Obs = (Obs[r,1,,,1]) / sum(Obs[r,1,,,1]) # Normalize observed temporary variable

      # Input into storage matrix
      Exp_mat[r,,] = array(tmp_Exp, dim = c(n_bins, n_sexes))
      Obs_mat[r,,] = array(tmp_Obs, dim = c(n_bins, n_sexes))
    } # end r loop
  } # end if 'Joint' comps by sex

  if(Comp_Type == 3) {
    tmp_Exp = aperm(Exp, perm = c(3,4,1,2,5)) # Reformat expected values so it's ordered by ages, sexes, and then regions
    tmp_Obs = aperm(Obs, perm = c(3,4,1,2,5)) # Reformat observed values so it's ordered by ages, sexes, and then regions

    # Expected values
    if(age_or_len == 0) { # if ages
      tmp_Exp = t(as.vector((tmp_Exp) / sum(tmp_Exp))) %*% kronecker(diag(n_regions * n_sexes), AgeingError) # apply ageing error
      tmp_Exp = as.vector((tmp_Exp)/ sum(tmp_Exp)) # renormalize to make sure sum to 1
    } # if ages

    if(age_or_len == 1) tmp_Exp = as.vector((Exp) / sum(Exp)) # Normalize temporary variable (lengths)

    tmp_Obs = (tmp_Obs) / sum(tmp_Obs) # Normalize observed temporary variable

    # Input into storage matrix
    Exp_mat[] = aperm(array(tmp_Exp, dim = c(n_bins, n_sexes, n_regions)), perm = c(3,1,2))
    Obs_mat[] = aperm(array(tmp_Obs, dim = c(n_bins, n_sexes, n_regions)), perm = c(3,1,2))
  } # Joint by region and sex

  return(list(Exp = Exp_mat, Obs = Obs_mat))

} # end function


#' Gets composition data proportions normalized according to the assessment specifications from RTMB
#'
#' @param data list of data inputs
#' @param rep report file from RTMB
#' @param year_labels vector of years
#' @param age_labels vector of age labels in assessment
#' @param len_labels vector of length labels in assessment
#'
#' @import dplyr
#' @importFrom tidyr drop_na
#' @returns List of fishery age, lengths, survey age, lengths dataframe as well as in matrix form (dimensioned by region, year, bin, sex, fleet)
#' @export get_comp_prop
#'
#' @examples
#' \dontrun{
#' comp_props <- get_comp_prop(data = data, rep = rep, age_labels = 2:31, len_labels = seq(41, 99, 2), year_labels = 1960:2024)
#' comp_props$Fishery_Ages %>%
#'   filter(Fleet == 1, Sex == 1) %>%
#'   ggplot() +
#'   geom_col(aes(x = Age, y = obs)) +
#'   geom_line(aes(x = Age, y = pred)) +
#'   facet_wrap(~Year, ncol = 3)
#'
#'   comp_props$Survey_Ages %>%
#'     group_by(Region, Age, Sex, Fleet) %>%
#'     summarize(lwr_obs = quantile(obs, 0.1),
#'               upr_obs = quantile(obs, 0.9),
#'               lwr_pred = quantile(pred, 0.1),
#'               upr_pred = quantile(pred, 0.9),
#'               obs = mean(obs),
#'               pred = mean(pred)) %>%
#'     ggplot() +
#'     geom_line(mapping = aes(x = Age, y = obs, color = 'Obs', lty = 'Obs'), lwd = 1.3) +
#'     geom_ribbon(mapping = aes(x = Age, y = obs, ymin = lwr_obs, ymax = upr_obs, fill = 'Obs'), alpha = 0.3) +
#'     geom_line(mapping = aes(x = Age, y = pred, color = 'Pred', lty = 'Pred'), lwd = 1.3) +
#'     geom_ribbon(mapping = aes(x = Age, y = pred, ymin = lwr_pred, ymax = upr_pred, fill = 'Pred'), alpha = 0.3) +
#'     facet_grid(Region~Fleet, labeller = labeller(
#'       Region = c('1' = "Region 1"),
#'       Fleet = c('1' = 'Domestic LL Survey', '3' = 'JP LL Survey')
#'     )) +
#'     labs(x = 'Age', y = 'Proportion', color = '', linetype = '', fill = '') +
#'     theme_bw(base_size = 20) +
#'     theme(legend.position = 'top')
#' }
get_comp_prop <- function(data,
                          rep,
                          age_labels,
                          len_labels,
                          year_labels
                          ) {

  # dimensinoing
  n_regions <- data$n_regions
  n_yrs <- length(data$years)
  n_ages <- length(data$ages)
  n_lens <- length(data$lens)
  n_sexes <- data$n_sexes
  n_fish_fleets <- data$n_fish_fleets
  n_srv_fleets <- data$n_srv_fleets

  # storage containers
  Obs_FishAge <- array(data = NA, dim = c(n_regions, n_yrs, n_ages, n_sexes, n_fish_fleets), dimnames = list(NULL, year_labels, age_labels, NULL, NULL)) # Obs fishery ages
  Obs_FishLen <- array(data = NA, dim = c(n_regions, n_yrs, n_lens, n_sexes, n_fish_fleets), dimnames = list(NULL, year_labels, len_labels, NULL, NULL)) # Obs fishery lengths
  Obs_SrvAge <- array(data = NA, dim = c(n_regions, n_yrs, n_ages, n_sexes, n_srv_fleets), dimnames = list(NULL, year_labels, age_labels, NULL, NULL)) # Obs survey ages
  Obs_SrvLen <- array(data = NA, dim = c(n_regions, n_yrs, n_lens, n_sexes, n_srv_fleets), dimnames = list(NULL, year_labels, len_labels, NULL, NULL)) # Obs survey lengths
  Pred_FishAge <- array(data = NA, dim = c(n_regions, n_yrs, n_ages, n_sexes, n_fish_fleets), dimnames = list(NULL, year_labels, age_labels, NULL, NULL)) # Predicted fishery ages
  Pred_FishLen <- array(data = NA, dim = c(n_regions, n_yrs, n_lens, n_sexes, n_fish_fleets), dimnames = list(NULL, year_labels, len_labels, NULL, NULL)) # Predicted fishery lengths
  Pred_SrvAge <- array(data = NA, dim = c(n_regions, n_yrs, n_ages, n_sexes, n_srv_fleets), dimnames = list(NULL, year_labels, age_labels, NULL, NULL)) # Predicted survey ages
  Pred_SrvLen <- array(data = NA, dim = c(n_regions, n_yrs, n_lens, n_sexes, n_srv_fleets), dimnames = list(NULL, year_labels, len_labels, NULL, NULL)) # Predicted survey lengths

  # Get quantities
  AgeingError <- data$AgeingError # ageing errors
  CAA <- rep$CAA # catch at age
  CAL <- rep$CAL # catch at len
  SrvIAA <- rep$SrvIAA # survey at age
  SrvIAL <- rep$SrvIAL # survey at length

  # Observed quantities
  ObsFishAgeComps <- data$ObsFishAgeComps
  ObsFishLenComps <- data$ObsFishLenComps
  ObsSrvAgeComps <- data$ObsSrvAgeComps
  ObsSrvLenComps <- data$ObsSrvLenComps

  # Composition Types
  FishAge_CompType <- data$FishAgeComps_Type
  SrvAge_CompType <- data$SrvAgeComps_Type
  FishLen_CompType <- data$FishLenComps_Type
  SrvLen_CompType <- data$SrvLenComps_Type

  # Aggregation Types
  FishAge_comp_agg_type <- data$FishAge_comp_agg_type
  FishLen_comp_agg_type <- data$FishLen_comp_agg_type
  SrvAge_comp_agg_type <- data$SrvAge_comp_agg_type
  SrvLen_comp_agg_type <- data$SrvLen_comp_agg_type

  # Whether ouse comp data
  UseFishAgeComps <- data$UseFishAgeComps
  UseFishLenComps <- data$UseFishLenComps
  UseSrvAgeComps <- data$UseSrvAgeComps
  UseSrvLenComps <- data$UseSrvLenComps

  # Fishery Ages
  for(y in 1:n_yrs) {
    for(f in 1:n_fish_fleets) {

      # figure out regions with obs
      use_regions <- which(UseFishAgeComps[,y,f] == 1)

      if(sum(use_regions) > 0) {
        Exp <- CAA[,y,,,f, drop = FALSE] # expected
        Obs <- ObsFishAgeComps[,y,,,f, drop = FALSE] # observed
        Comp_Type <- FishAge_CompType[y,f] # composition type
        comp_agg_type <- FishAge_comp_agg_type[f] # aggregation type

        # reformat expected compositions
        tmp_comps <- Restrc_Comps(Exp = Exp, Obs = Obs, Comp_Type = Comp_Type,
                                  age_or_len = 0, AgeingError = AgeingError, comp_agg_type = comp_agg_type)
        # Input into storage
        Obs_FishAge[,y,,,f] <- tmp_comps$Obs
        Pred_FishAge[,y,,,f] <- tmp_comps$Exp
      }

    } # end f
  } # end y

  # Fishery Lengths
  for(y in 1:n_yrs) {
    for(f in 1:n_fish_fleets) {

      use_regions <- which(UseFishLenComps[,y,f] == 1) # figure out regions with obs

      if(sum(use_regions) > 0) {
        Exp <- CAL[,y,,,f, drop = FALSE] # expected
        Obs <- ObsFishLenComps[,y,,,f, drop = FALSE] # observed
        Comp_Type <- FishLen_CompType[y,f] # composition type
        comp_agg_type <- FishLen_comp_agg_type[f] # aggregation type

        # get compositions
        tmp_comps <- Restrc_Comps(Exp = Exp, Obs = Obs, Comp_Type = Comp_Type,
                                  age_or_len = 1, AgeingError = NA, comp_agg_type = comp_agg_type)
        # Input into storage
        Obs_FishLen[,y,,,f] <- tmp_comps$Obs
        Pred_FishLen[,y,,,f] <- tmp_comps$Exp
      }

    } # end f
  } # end y

  # Survey Ages
  for(y in 1:n_yrs) {
    for(f in 1:n_srv_fleets) {

      # figure out regions with obs
      use_regions <- which(UseSrvAgeComps[,y,f] == 1)

      if(sum(use_regions) > 0) {
        Exp <- SrvIAA[,y,,,f, drop = FALSE] # expected
        Obs <- ObsSrvAgeComps[,y,,,f, drop = FALSE] # observed
        Comp_Type <- SrvAge_CompType[y,f] # composition type
        comp_agg_type <- SrvAge_comp_agg_type[f] # aggregation type

        # reformat expected compositions
        tmp_comps <- Restrc_Comps(Exp = Exp, Obs = Obs, Comp_Type = Comp_Type,
                                  age_or_len = 0, AgeingError = AgeingError, comp_agg_type = comp_agg_type)
        # Input into storage
        Obs_SrvAge[,y,,,f] <- tmp_comps$Obs
        Pred_SrvAge[,y,,,f] <- tmp_comps$Exp
      }

    } # end f
  } # end y

  # Survey Lengths
  for(y in 1:n_yrs) {
    for(f in 1:n_srv_fleets) {

      # figure out regions with obs
      use_regions <- which(UseSrvLenComps[,y,f] == 1)

      if(sum(use_regions) > 0) {
        Exp <- SrvIAL[,y,,,f, drop = FALSE] # expected
        Obs <- ObsSrvLenComps[,y,,,f, drop = FALSE] # observed
        Comp_Type <- SrvLen_CompType[y,f] # composition type
        comp_agg_type <- SrvLen_comp_agg_type[f] # aggregation type

        # reformat expected compositions
        tmp_comps <- Restrc_Comps(Exp = Exp, Obs = Obs, Comp_Type = Comp_Type,
                                  age_or_len = 1, AgeingError = NA, comp_agg_type = comp_agg_type)
        # Input into storage
        Obs_SrvLen[,y,,,f] <- tmp_comps$Obs
        Pred_SrvLen[,y,,,f] <- tmp_comps$Exp
      }

    } # end f
  } # end y

  # Process outputs
  all_fishages <- reshape2::melt(Obs_FishAge) %>%
    dplyr::rename(obs = value) %>%
    tidyr::drop_na() %>%
    dplyr::left_join(reshape2::melt(Pred_FishAge) %>% dplyr::rename(pred = value), by = c("Var1", "Var2", "Var3", "Var4", "Var5")) %>%
    dplyr::rename(Region = Var1, Year = Var2, Age = Var3, Sex = Var4, Fleet = Var5) %>%
    dplyr::mutate(Type = 'Fishery Ages')

  # Fish Lengths
  all_fishlens <- reshape2::melt(Obs_FishLen) %>%
    dplyr::rename(obs = value) %>%
    tidyr::drop_na() %>%
    dplyr::left_join(reshape2::melt(Pred_FishLen) %>% dplyr::rename(pred = value), by = c("Var1", "Var2", "Var3", "Var4", "Var5")) %>%
    dplyr::rename(Region = Var1, Year = Var2, Len = Var3, Sex = Var4, Fleet = Var5) %>%
    dplyr::mutate(Type = 'Fishery Lengths')

  # Survey Ages
  all_srvages <- reshape2::melt(Obs_SrvAge) %>%
    dplyr::rename(obs = value) %>%
    tidyr::drop_na() %>%
    dplyr::left_join(reshape2::melt(Pred_SrvAge) %>% dplyr::rename(pred = value), by = c("Var1", "Var2", "Var3", "Var4", "Var5")) %>%
    dplyr::rename(Region = Var1, Year = Var2, Age = Var3, Sex = Var4, Fleet = Var5) %>%
    dplyr::mutate(Type = 'Survey Ages')

  # Survey Lengths
  all_srvlens <- reshape2::melt(Obs_SrvLen) %>%
    dplyr::rename(obs = value) %>%
    tidyr::drop_na() %>%
    dplyr::left_join(reshape2::melt(Pred_SrvLen) %>% dplyr::rename(pred = value), by = c("Var1", "Var2", "Var3", "Var4", "Var5")) %>%
    dplyr::rename(Region = Var1, Year = Var2, Len = Var3, Sex = Var4, Fleet = Var5) %>%
    dplyr::mutate(Type = 'Survey Lengths')

  return(list(# data frames of observed and expected comps
              Fishery_Ages = all_fishages,
              Fishery_Lens = all_fishlens,
              Survey_Ages = all_srvages,
              Survey_Lens = all_srvlens,

              # Arrays of observed comps
              Obs_FishAge_mat = Obs_FishAge,
              Obs_FishLen_mat = Obs_FishLen,
              Obs_SrvAge_mat = Obs_SrvAge,
              Obs_SrvLen_mat = Obs_SrvLen,

              # Arrays of expected comps
              Pred_FishAge_mat = Pred_FishAge,
              Pred_FishLen_mat = Pred_FishLen,
              Pred_SrvAge_mat = Pred_SrvAge,
              Pred_SrvLen_mat = Pred_SrvLen))
} # end function


#' Function to format comp data and get OSA residuals from afscOSA (uses afscOSA and compresid as backend to get OSA residuals)
#'
#' @param obs_mat Matrix of observed values, which can have NAs - gets removed with years arge (dimensioned by region, year, age, sex, fleet)
#' @param exp_mat Matrix of expceted values, which can have NAs - gets removed with years arg (dimensioned by region, year, age, sex, fleet)
#' @param N Input or effective sample size
#' @param years Years we want to point to and filter to
#' @param fleet Fleet we want to filter to
#' @param bins Vector of age or length bins
#' @param comp_type Composition type - whether this is aggregated == 0, split by region and sex == 1, split by region joint by sex = =2, and joint by region and sex == 3
#' @param bin_label Bin label for whether these are ages or lengths
#'
#' @import dplyr
#' @import compResidual
#' @import afscOSA
#' @returns Dataframe of OSA residuals
#' @export get_osa
#'
#' @examples
#' \dontrun{
#' comp_props <- get_comp_prop(data = data, rep = rep, age_labels = 2:31, len_labels = seq(41, 99, 2), year_labels = 1960:2024)
#' osa_results <- get_osa(obs_mat = comp_props$Obs_FishAge_mat,
#'                       exp_mat = comp_props$Pred_FishAge_mat,
#'                       N = 20 * data$Wt_FishAgeComps[1,1,1],
#'                       years = 1999:2023,
#'                       fleet = 1,
#'                       bins = 2:31,
#'                       comp_type = 0,
#'                       bin_label = "Age")
#' }
get_osa <- function(obs_mat,
                    exp_mat,
                    N,
                    years,
                    fleet,
                    bins,
                    comp_type,
                    bin_label
                    ) {

  years <- as.character(years) # define as character
  obs <- obs_mat[,years,,,fleet, drop = FALSE] # get filtered observed matrix
  exp <- exp_mat[,years,,,fleet, drop = FALSE] # get filtered expected matrix
  n_regions <- dim(obs)[1]
  n_sexes <- dim(obs)[4]

  # if comps are aggregated
  if(comp_type == 0) {
    tmp_obs <- obs[1,,,1,1] # only get a single sex out
    tmp_exp <- exp[1,,,1,1] # only get a single sex out

    # compute OSA
    tmp_osa <- afscOSA::run_osa(obs = tmp_obs, exp = tmp_exp, N = N, years = years,
                   index = bins, fleet = as.character(fleet), index_label = bin_label)

    tmp_osa$res$comp_type <- "Aggregated"
    tmp_osa$agg$comp_type <- "Aggregated"
    osa_all <- tmp_osa
  }

  # if comps are split by region and sex
  if(comp_type == 1) {

    # empty dataframes to bind to
    res_all <- data.frame()
    agg_all <- data.frame()

    for(r in 1:n_regions) {
      for(s in 1:n_sexes) {
        tmp_obs <- obs[r,,,s,1] # get observations
        tmp_exp <- exp[r,,,s,1] # get expected

        # compute OSA
        tmp_osa <- afscOSA::run_osa(obs = tmp_obs, exp = tmp_exp, N = N, years = years,
                       index = bins, fleet = as.character(fleet), index_label = bin_label)

        # Doing some naming stuff
        tmp_osa$res$region <- r
        tmp_osa$agg$region <- r
        tmp_osa$res$sex <- s
        tmp_osa$agg$sex <- s
        tmp_osa$res$comp_type <- "SpltR_SpltS"
        tmp_osa$agg$comp_type <- "SpltR_SpltS"

        res_all <- rbind(res_all, tmp_osa$res)
        agg_all <- rbind(agg_all, tmp_osa$agg)

      } # end s loop
    } # end r loop

    osa_all <- list(res = res_all, agg = agg_all)

  } # end split region and sex

  # if comp types are split by sex, joint by region
  if(comp_type == 2) {

    # empty dataframes to bind to
    res_all <- data.frame()
    agg_all <- data.frame()

    for(r in 1:n_regions) {

      # initialize to cbind
      tmp_obs <- NULL
      tmp_exp <- NULL

      for(s in 1:n_sexes) {
        tmp_obs <- cbind(tmp_obs, obs[r,,,s,1]) # get observations
        tmp_exp <- cbind(tmp_exp, exp[r,,,s,1]) # get expected
      } # end s loop

      # compute OSA
      tmp_osa <- afscOSA::run_osa(obs = tmp_obs, exp = tmp_exp, N = N, years = years,
                         index = paste(rep(1:n_sexes, each = length(bins)), "_", rep(bins, times = n_sexes), sep = ""),
                         fleet = as.character(fleet), index_label = bin_label)

        # Doing some naming stuff
        tmp_osa$res$region <- r
        tmp_osa$agg$region <- r
        tmp_osa$res <- tmp_osa$res %>% dplyr::mutate(split_index = str_split(index, "_"),  # Split once and store as list
                                              sex = sapply(split_index, `[`, 1),
                                              index = sapply(split_index, `[`, 2)) %>% dplyr::select(-split_index)
        tmp_osa$agg <- tmp_osa$agg %>% dplyr::mutate(split_index = str_split(index, "_"),  # Split once and store as list
                                             sex = sapply(split_index, `[`, 1),
                                             index = sapply(split_index, `[`, 2)) %>% dplyr::select(-split_index)
        tmp_osa$res$comp_type <- "SpltR_JntS"
        tmp_osa$agg$comp_type <- "SpltR_JntS"

        res_all <- rbind(res_all, tmp_osa$res)
        agg_all <- rbind(agg_all, tmp_osa$agg)

    } # end r loop

    osa_all <- list(res = res_all, agg = agg_all)

  } # end split region, joint by sex

  if(comp_type == 3) {

    # initialize to cbind
    tmp_obs <- NULL
    tmp_exp <- NULL

    for(r in 1:n_regions) {
      for(s in 1:n_sexes) {
        tmp_obs <- cbind(tmp_obs, obs[r,,,s,1]) # get observations
        tmp_exp <- cbind(tmp_exp, exp[r,,,s,1]) # get expected
      } # end s loop
    } # end r loop

    # compute OSA
    tmp_osa <- afscOSA::run_osa(obs = tmp_obs, exp = tmp_exp, N = N, years = years,
                       index = paste(rep(1:n_sexes, each = length(bins) * n_regions), "_",
                                     rep(bins, each = n_regions, times = n_sexes), "_",
                                     rep(1:n_regions, times = length(bins) * n_sexes),
                                     sep = ""),
                       fleet = as.character(fleet), index_label = bin_label)

    # Doing some naming stuff
    tmp_osa$res <- tmp_osa$res %>%
      dplyr::mutate(split_index = str_split(index, "_"),  # Split once and store as list
              sex = sapply(split_index, `[`, 1),
              index = sapply(split_index, `[`, 2),
              region = sapply(split_index, `[`, 3)) %>%
      dplyr::select(-split_index)

    tmp_osa$agg <- tmp_osa$agg %>%
      dplyr::mutate(split_index = str_split(index, "_"),  # Split once and store as list
             sex = sapply(split_index, `[`, 1), # get sex
             index = sapply(split_index, `[`, 2), # get bin
             region = sapply(split_index, `[`, 3)) %>% # get region
      dplyr::select(-split_index)

    tmp_osa$res$comp_type <- "JntR_JntS"
    tmp_osa$agg$comp_type <- "JntR_JntS"

    osa_all <- tmp_osa

  } # end joint region and sex

  return(osa_all)
}

#' Plots OSA residuals from outputs from get_osa. Much of this code is taken from the afscOM package, but with modificaitons to plot features.
#'
#' @param osa_results List object obtained from get_osa, that contains a dataframe of residuals and aggregated fits.
#'
#' @returns A vareity of plots for OSA residuals (list)
#' @export plot_resids
#' @import dplyr
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#' comp_props <- get_comp_prop(data = data, rep = rep, age_labels = 2:31, len_labels = seq(41, 99, 2), year_labels = 1960:2024)
#' osa_results <- get_osa(obs_mat = comp_props$Obs_FishLen_mat,
#'                        exp_mat = comp_props$Pred_FishLen_mat,
#'                        N = 20 * data$Wt_FishAgeComps[1,1,1],
#'                        years = 1999:2023,
#'                        fleet = 1,
#'                        bins = 2:31,
#'                        comp_type = 1,
#'                        bin_label = "Age")
#'
#' osa_plot <- plot_resids(osa_results)
#' }
plot_resids <- function(osa_results) {

  # extract results
  res <- osa_results$res %>% dplyr::mutate(sign = ifelse(resid < 0, "Neg", "Pos"), Outlier = ifelse(abs(resid) > 3, "Yes", "No"))
  agg <- osa_results$agg

  # Aggregated Comps
  if(unique(res$comp_type) == "Aggregated") {

    # Get standarized normal residuals
    sdnr <- res %>% dplyr::summarise(sdnr = paste0("SDNR = ", formatC(round(sd(resid), 3), format = "f", digits = 2)))

    # sdnr plot
    sdnr_plot <- ggplot() +
      geom_abline(slope = 1, intercept = 0, lty = 2, lwd = 1.3) +
      stat_qq(data = res, aes(sample = resid), col = "blue", size = 2, alpha = 0.5) +
      theme_bw(base_size = 20) +
      geom_text(data = sdnr, aes(x = -Inf, y = Inf, label = sdnr), hjust = -0.5, vjust = 2.5, size = 8)
  }

  # Split Sex and Split Region
  if(unique(res$comp_type) == "SpltR_SpltS") {

    # Get standarized normal residuals
    sdnr <- res %>% dplyr::group_by(region, sex) %>%
      dplyr::summarise(sdnr = paste0("SDNR = ", formatC(round(sd(resid), 3), format = "f", digits = 2)))

    # sdnr plot
    sdnr_plot <- ggplot() +
      geom_abline(slope = 1, intercept = 0, lty = 2, lwd = 1.3) +
      stat_qq(data = res, aes(sample = resid), col = "blue", size = 2, alpha = 0.5) +
      labs(x = "Theoretical quantiles", y = "Sample quantiles") +
      facet_grid(region ~ sex, labeller = labeller(
        region = function(x) paste0("Region ", x),
        sex = function(x) paste0("Sex ", x)
      )) +
      theme_bw(base_size = 20) +
      geom_text(data = sdnr, aes(x = -Inf, y = Inf, label = sdnr), hjust = -0.5, vjust = 2.5, size = 8)
  }


  # Joint Sex and Split Region
  if(unique(res$comp_type) == "SpltR_JntS") {

    # Get standarized normal residuals
    sdnr <- res %>% dplyr::group_by(region) %>%
      dplyr::summarise(sdnr = paste0("SDNR = ", formatC(round(sd(resid), 3), format = "f", digits = 2)))

    # sdnr plot
    sdnr_plot <- ggplot() +
      geom_abline(slope = 1, intercept = 0, lty = 2, lwd = 1.3) +
      stat_qq(data = res, aes(sample = resid), col = "blue", size = 2, alpha = 0.5) +
      labs(x = "Theoretical quantiles", y = "Sample quantiles") +
      facet_grid(region ~ sex, labeller = labeller(
        region = function(x) paste0("Region ", x),
        sex = function(x) paste0("Sex ", x)
      )) +
      theme_bw(base_size = 20) +
      geom_text(data = sdnr, aes(x = -Inf, y = Inf, label = sdnr), hjust = -0.5, vjust = 2.5, size = 8)
  }

  # Joint Sex and Split Region
  if(unique(res$comp_type) == "JntR_JntS") {

    # Get standarized normal residuals
    sdnr <- res %>% dplyr::summarise(sdnr = paste0("SDNR = ", formatC(round(sd(resid), 3), format = "f", digits = 2)))

    # sdnr plot
    sdnr_plot <- ggplot() +
      geom_abline(slope = 1, intercept = 0, lty = 2, lwd = 1.3) +
      stat_qq(data = res, aes(sample = resid), col = "blue", size = 2, alpha = 0.5) +
      labs(x = "Theoretical quantiles", y = "Sample quantiles") +
      facet_grid(region ~ sex, labeller = labeller(
        region = function(x) paste0("Region ", x),
        sex = function(x) paste0("Sex ", x)
      )) +
      theme_bw(base_size = 20) +
      geom_text(data = sdnr, aes(x = -Inf, y = Inf, label = sdnr), hjust = -0.5, vjust = 2.5, size = 8)
  }

  # bubble plot
  bubble_plot <- ggplot(data = res, aes(x = year, y = as.numeric(index),
                                        color = sign, size = abs(resid), shape = Outlier, alpha = abs(resid))) +
    geom_point() +
    scale_color_manual(values = c("blue", "red")) +
    labs(x = "Year", y = unique(res$index_label), color = "Sign",
         sign = "abs(Resid)", size = "abs(Resid)", alpha = "abs(Resid)") +
    facet_grid(region ~ sex, labeller = labeller(
      region = function(x) paste0("Region ", x),
      sex = function(x) paste0("Sex ", x)
    )) +
    theme_bw(base_size = 20) +
    theme(legend.position = 'top')

  # aggregated plot
  agg_plot <- ggplot(data = agg) +
    geom_bar(aes(x = as.numeric(index), y = obs), stat = "identity", color = "blue", fill = "blue", alpha = 0.4) +
    geom_point(aes(x = as.numeric(index), y = exp), color = "red") +
    geom_line(aes(x = as.numeric(index), y = exp), color = "red") +
    facet_grid(region ~ sex, labeller = labeller(
      region = function(x) paste0("Region ", x),
      sex = function(x) paste0("Sex ", x)
    )) +
    theme_bw(base_size = 20) +
    theme(legend.position = 'top') +
    labs(x = unique(res$index_label), y = "Proportion")

  return(list(sdnr_plot, bubble_plot, agg_plot))

}

