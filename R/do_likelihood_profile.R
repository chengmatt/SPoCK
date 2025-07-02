#' Run Likelihood Profile
#'
#' @param data data list from model
#' @param parameters parameter list from model
#' @param mapping mapping list from model
#' @param random character vector of random effects to estimate
#' @param what parameter name we want to profile
#' @param idx Index for an parameter array, pointing to the value we want to map off (index is relative to a flattened array)
#' @param min_val minimum value of profile
#' @param max_val maximum value of profile
#' @param inc increment value between min and max value
#'
#' @import dplyr
#' @import RTMB
#' @importFrom reshape2 melt
#' @importFrom stats rnorm nlminb
#' @returns Returns a list of likelihood profiled values for each data component with their respective dimensions (e.g., likelihood profiles by fleet, region, year, etc.) as well likelihood profiles for each data component, aggregated across all their respective dimensions.
#' @export do_likelihood_profile
#'
do_likelihood_profile <- function(data,
                                  parameters,
                                  mapping,
                                  random = NULL,
                                  what,
                                  idx = NULL,
                                  min_val,
                                  max_val,
                                  inc = 0.05
                                  ) {

  if(min_val > max_val) {
    stop("`min_val` is greater than `max_val`. This likely occurred because you are profiling a log-transformed parameter. Try swapping the values: use the current `min_val` as `max_val`, and vice versa.")
  }

  # create values to profile across
  vals <- seq(min_val, max_val, inc)

  # Create objects to store values
  jnLL <- matrix(NA, nrow = length(vals), ncol = 1, dimnames = list(vals, NULL))
  rec_nLL <- matrix(NA, nrow = length(vals), ncol = 1, dimnames = list(vals, NULL))
  sel_nLL <- matrix(NA, nrow = length(vals), ncol = 1, dimnames = list(vals, NULL))
  M_nLL <- matrix(NA, nrow = length(vals), ncol = 1, dimnames = list(vals, NULL))
  h_nLL <- matrix(NA, nrow = length(vals), ncol = 1, dimnames = list(vals, NULL))
  Movement_nLL <- matrix(NA, nrow = length(vals), ncol = 1, dimnames = list(vals, NULL))
  TagRep_nLL <- matrix(NA, nrow = length(vals), ncol = 1, dimnames = list(vals, NULL))
  Fmort_nLL <- matrix(NA, nrow = length(vals), ncol = 1, dimnames = list(vals, NULL))
  Tag_nLL <- data.frame()
  Catch_nLL <- data.frame()
  FishAge_nLL <- data.frame()
  SrvAge_nLL <- data.frame()
  SrvLen_nLL <- data.frame()
  FishLen_nLL <- data.frame()
  FishIdx_nLL <- data.frame()
  SrvIdx_nLL <- data.frame()

  # If there is more than one value in this parameter
  for(j in 1:length(vals)) {
    if(!is.null(dim(parameters[[what]]))) {

      # Input fixed value into parameter list
      parameters[[what]] <- do.call(`[<-`, c(list(parameters[[what]]), idx, list(vals[j])))

      # Now, figure out which parameter to map off
      counter <- 1
      map_parameter <- do.call(`[<-`, c(list(parameters[[what]]), idx, list(NA))) # input NA
      non_na <- which(!is.na(map_parameter)) # figure out non NA positions and fill with unique numbers
      for(i in 1:length(non_na)) {
        map_parameter[non_na[i]] <- counter
        counter <- counter + 1 # update counter
      } # end i
      mapping[[what]] <- factor(map_parameter) # input NA into mapping list

    } else { # else, there is only one value in this parameter
      parameters[[what]] <- vals[j]
      mapping[[what]] <- factor(NA)
    }

    # make adfun
    SPoCK_rtmb_model <- RTMB::MakeADFun(cmb(SPoCK_rtmb, data), parameters = parameters, map = mapping, random = random, silent = TRUE)

    # Within your loop
    tryCatch({
      SPoCK_optim <- stats::nlminb(SPoCK_rtmb_model$par, SPoCK_rtmb_model$fn, SPoCK_rtmb_model$gr,
                                   control = list(iter.max = 1e6, eval.max = 1e6, rel.tol = 1e-15))

      # Get report
      report <- SPoCK_rtmb_model$report(SPoCK_rtmb_model$env$last.par.best)

      # Store values and save
      jnLL[j,1] <- report$jnLL
      rec_nLL[j,1] <- sum(report$Init_Rec_nLL, report$Rec_nLL)
      M_nLL[j,1] <- report$M_nLL
      sel_nLL[j,1] <- report$sel_nLL
      Movement_nLL[j,1] <- report$Movement_nLL
      h_nLL[j,1] <- report$h_nLL
      TagRep_nLL[j,1] <- report$TagRep_nLL
      Fmort_nLL[j,1] <- sum(report$Fmor_nLL)
      Tag_nLL <- rbind(Tag_nLL, reshape2::melt(report$Tag_nLL) %>% dplyr::mutate(prof_val = vals[j]))
      Catch_nLL <- rbind(Catch_nLL, reshape2::melt(report$Catch_nLL) %>% dplyr::mutate(prof_val = vals[j]))
      FishAge_nLL <- rbind(FishAge_nLL, reshape2::melt(report$FishAgeComps_nLL) %>% dplyr::mutate(prof_val = vals[j]))
      SrvAge_nLL <- rbind(SrvAge_nLL, reshape2::melt(report$SrvAgeComps_nLL) %>% dplyr::mutate(prof_val = vals[j]))
      SrvLen_nLL <- rbind(SrvLen_nLL, reshape2::melt(report$SrvLenComps_nLL) %>% dplyr::mutate(prof_val = vals[j]))
      FishLen_nLL <- rbind(FishLen_nLL, reshape2::melt(report$FishLenComps_nLL) %>% dplyr::mutate(prof_val = vals[j]))
      FishIdx_nLL <- rbind(FishIdx_nLL, reshape2::melt(report$FishIdx_nLL) %>% dplyr::mutate(prof_val = vals[j]))
      SrvIdx_nLL <- rbind(SrvIdx_nLL, reshape2::melt(report$SrvIdx_nLL) %>% dplyr::mutate(prof_val = vals[j]))

      print(paste("Likelihood profile is at:", j / length(vals) * 100))

    }, error = function(e) {
      message("Failed to optimize: ", e$message)
    })

  } # end j loop

  # Doing some residual munging into the correct format
  jnLL_df <- reshape2::melt(jnLL) %>%
    dplyr::select(-Var2) %>%
    dplyr::rename(prof_val = Var1) %>%
    dplyr::mutate(type = 'jnLL')
  rec_nLL_df <- reshape2::melt(rec_nLL) %>%
    dplyr::select(-Var2) %>%
    dplyr::rename(prof_val = Var1) %>%
    dplyr::mutate(type = 'RecPen')
  M_nLL_df <- reshape2::melt(M_nLL) %>%
    dplyr::select(-Var2) %>%
    dplyr::rename(prof_val = Var1) %>%
    dplyr::mutate(type = 'M Prior')
  sel_nLL_df <- reshape2::melt(sel_nLL) %>%
    dplyr::select(-Var2) %>%
    dplyr::rename(prof_val = Var1) %>%
    dplyr::mutate(type = 'Selex Pen')
  Movement_nLL_df <- reshape2::melt(Movement_nLL) %>%
    dplyr::select(-Var2) %>%
    dplyr::rename(prof_val = Var1) %>%
    dplyr::mutate(type = 'Move Prior')
  h_nLL_df <- reshape2::melt(h_nLL) %>% dplyr::select(-Var2) %>% dplyr::rename(prof_val = Var1) %>% dplyr::mutate(type = 'h Prior')
  TagRep_nLL_df <- reshape2::melt(TagRep_nLL) %>%
    dplyr::select(-Var2) %>%
    dplyr::rename(prof_val = Var1) %>%
    dplyr::mutate(type = 'TagRep Prior')
  Fmort_nLL_df <- reshape2::melt(Fmort_nLL) %>%
    dplyr::select(-Var2) %>%
    dplyr::rename(prof_val = Var1) %>%
    dplyr::mutate(type = 'FmortPen')
  Catch_nLL_df <- Catch_nLL %>%
    dplyr::rename(Region = Var1, Year = Var2, Fleet = Var3) %>%
    dplyr::mutate(type = 'Catch')
  FishAge_nLL_df <- FishAge_nLL %>%
    dplyr::rename(Region = Var1, Year = Var2, Sex = Var3, Fleet = Var4) %>%
    dplyr::mutate(type = 'FishAge')
  SrvAge_nLL_df <- SrvAge_nLL %>%
    dplyr::rename(Region = Var1, Year = Var2, Sex = Var3, Fleet = Var4) %>%
    dplyr::mutate(type = 'SrvAge')
  FishLen_nLL_df <- FishLen_nLL %>%
    dplyr::rename(Region = Var1, Year = Var2, Sex = Var3, Fleet = Var4) %>%
    dplyr::mutate(type = 'FishLen')
  SrvLen_nLL_df <- SrvLen_nLL %>%
    dplyr::rename(Region = Var1, Year = Var2, Sex = Var3, Fleet = Var4) %>%
    dplyr::mutate(type = 'SrvLen')
  FishIdx_nLL_df <- FishLen_nLL %>%
    dplyr::rename(Region = Var1, Year = Var2, Fleet = Var3) %>%
    dplyr::mutate(type = 'FishIdx')
  SrvIdx_nLL_df <- SrvIdx_nLL %>%
    dplyr::rename(Region = Var1, Year = Var2, Fleet = Var3) %>%
    dplyr::mutate(type = 'SrvIdx')
  Tag_nLL_df <- Tag_nLL %>%
    dplyr::rename(Recap_Year = Var1, Tag_Cohort = Var2, Region = Var3, Age = Var4, Sex = Var5) %>%
    dplyr::mutate(type = 'Tagging')

  # Get likelihoods aggregated across all dimensions
  agg_nLL <- rbind(jnLL_df, rec_nLL_df, M_nLL_df, Movement_nLL_df, h_nLL_df,
                   TagRep_nLL_df,Fmort_nLL_df, sel_nLL_df,
                   Catch_nLL_df %>% dplyr::group_by(prof_val, type) %>%
                     dplyr::summarize(value = sum(value)),
                   Tag_nLL_df %>% dplyr::group_by(prof_val, type) %>%
                     dplyr::summarize(value = sum(value, na.rm = T)),
                   FishAge_nLL_df %>% dplyr::group_by(prof_val, type) %>%
                     dplyr::summarize(value = sum(value, na.rm = T)),
                   SrvAge_nLL_df %>% dplyr::group_by(prof_val, type) %>%
                     dplyr::summarize(value = sum(value, na.rm = T)),
                   FishLen_nLL_df %>% dplyr::group_by(prof_val, type) %>%
                     dplyr::summarize(value = sum(value, na.rm = T)),
                   SrvLen_nLL_df %>% dplyr::group_by(prof_val, type) %>%
                     dplyr::summarize(value = sum(value, na.rm = T)),
                   FishIdx_nLL_df %>% dplyr::group_by(prof_val, type) %>%
                     dplyr::summarize(value = sum(value, na.rm = T)),
                   SrvIdx_nLL_df %>% dplyr::group_by(prof_val, type) %>%
                     dplyr::summarize(value = sum(value, na.rm = T))
                   )

  profile_list <- list(jnLL_df = jnLL_df,
                       rec_nLL_df = rec_nLL_df,
                       M_nLL_df = M_nLL_df,
                       sel_nLL_df = sel_nLL_df,
                       Movement_nLL_df = Movement_nLL_df,
                       h_nLL_df = h_nLL_df,
                       TagRep_nLL_df = TagRep_nLL_df,
                       Fmort_nLL_df = Fmort_nLL_df,
                       Catch_nLL_df = Catch_nLL_df,
                       Tag_nLL_df = Tag_nLL_df,
                       FishAge_nLL_df = FishAge_nLL_df,
                       SrvAge_nLL_df = SrvAge_nLL_df,
                       FishLen_nLL_df = FishLen_nLL_df,
                       SrvLen_nLL_df = SrvLen_nLL_df,
                       FishIdx_nLL_df = FishIdx_nLL_df,
                       SrvIdx_nLL_df = SrvIdx_nLL_df,
                       agg_nLL = agg_nLL
                       )

  return(profile_list)
}





