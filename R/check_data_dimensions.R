#' Helper function used to check data dimensions to ensure they are correct
#'
#' @param x Object to evaluate
#' @param n_regions Number of regions
#' @param n_years Number of years
#' @param n_ages Number of ages
#' @param n_lens Number of lengths
#' @param n_sexes Number of sexes
#' @param n_fish_fleets Number of fishery fleets
#' @param n_srv_fleets Number of survey fleets
#' @param what charcter specifying what to be evaluated
#'
#' @export check_data_dimensions
#'
check_data_dimensions <- function(x,
                                  n_regions = NULL,
                                  n_years = NULL,
                                  n_ages = NULL,
                                  n_lens = NULL,
                                  n_sexes = NULL,
                                  n_fish_fleets = NULL,
                                  n_srv_fleets = NULL,
                                  what
                                  ) {

# Biologicals -------------------------------------------------------------

  # Weight at age (spawning), maturity, or natural mortality
  if(what %in% c('WAA', 'MatAA', "Fixed_natmort")) {
    if(sum(dim(x) == c(n_regions, n_years, n_ages, n_sexes)) != 4)
      stop(paste("Dimensions of", what, "are not correct. Should be n_regions, n_years, n_ages, and n_sexes"))
  }

  # weight at age for the fishery
  if(what == 'WAA_fish') {
    if(sum(dim(x) == c(n_regions, n_years, n_ages, n_sexes, n_fish_fleets)) != 5)
      stop(paste("Dimensions of", what, "are not correct. Should be n_regions, n_years, n_ages, n_sexes, and n_fish_fleets"))
  }

  # weight at age for the survey
  if(what == 'WAA_srv') {
    if(sum(dim(x) == c(n_regions, n_years, n_ages, n_sexes, n_srv_fleets)) != 5)
      stop(paste("Dimensions of", what, "are not correct. Should be n_regions, n_years, n_ages, n_sexes, and n_srv_fleets"))
  }

  if(what == 'AgeingError') { # Not checking the age dimension
    if(sum(dim(x)[1] == n_ages) != 1)
      stop("Dimensions of AgeingError are not correct. Should be n_ages, number of observed composition ages")
  }

  if(what == 'SizeAgeTrans') {
    if(sum(dim(x) == c(n_regions, n_years, n_lens, n_ages, n_sexes)) != 5)
       stop("Dimensions of SizeAgeTrans are not correct. Should be n_regions, n_years, n_lens, n_ages, and n_sexes")
  }

  if(what == 'Fixed_Movement') {
    if(sum(dim(x) == c(n_regions, n_regions, n_years, n_ages, n_sexes)) != 5)
      stop("Fixed Movement Matrix does not have the correct dimensions. This should be n_regions, n_regions, n_years, n_ages, n_sexes")
  }


# Fishery Stuff -----------------------------------------------------------------

  if(what %in% c('ObsCatch', "UseCatch", 'ObsFishIdx', 'ObsFishIdx_SE', 'UseFishIdx', 'UseFishAgeComps', 'UseFishLenComps')) {
    if(sum(dim(x) == c(n_regions, n_years, n_fish_fleets)) != 3)
      stop(paste(what, " is not the correct dimension. Should be n_regions, n_years, n_fish_fleets"))
  }

  if(what == 'Catch_Type') {
    if(sum(dim(x) == c(n_years, n_fish_fleets)) != 2)
      stop("Catch_Type is not the correct dimension. Should be n_years, n_fish_fleets")
  }

  if(what %in% c('Catch_Constant', "fish_idx_type", "FishAgeComps_LikeType", "FishLenComps_LikeType")) {
    if(length(x) != n_fish_fleets)
      stop(paste(what, "needs to have a length of n_fish_fleets"))
  }

  if(what == "ObsFishAgeComps") { # Not checking the age dimension
    if(sum(dim(x)[-3] == c(n_regions, n_years, n_sexes, n_fish_fleets)) != 4)
      stop(paste("ObsFishAgeComps is not the correct dimension. Should be n_regions, n_years, number of observed composition ages, n_sexes, n_fish_fleets"))
  }

  if(what == "ObsFishLenComps") {
    if(sum(dim(x) == c(n_regions, n_years, n_lens, n_sexes, n_fish_fleets)) != 5)
      stop(paste("ObsFishLenComps is not the correct dimension. Should be n_regions, n_years, n_lens, n_sexes, n_fish_fleets"))
  }

  if(what %in% c('ISS_FishLenComps', 'ISS_FishAgeComps')) {
    if(sum(dim(x) == c(n_regions, n_years, n_sexes, n_fish_fleets)) != 4)
      stop(paste(what, " is not the correct dimension. Should be n_regions, n_years, n_sexes, n_fish_fleets"))
  }


# Survey Stuff ------------------------------------------------------------

  if(what %in% c('ObsSrvIdx', 'ObsSrvIdx_SE', 'UseSrvIdx', 'UseSrvAgeComps', 'UseSrvLenComps')) {
    if(sum(dim(x) == c(n_regions, n_years, n_srv_fleets)) != 3)
      stop(paste(what, " is not the correct dimension. Should be n_regions, n_years, n_srv_fleets"))
  }

  if(what %in% c("srv_idx_type", "SrvAgeComps_LikeType", "SrvLenComps_LikeType")) {
    if(length(x) != n_srv_fleets)
      stop(paste(what, "needs to have a length of n_srv_fleets"))
  }

  if(what == "ObsSrvAgeComps") { # Not checking the age dimension
    if(sum(dim(x)[-3] == c(n_regions, n_years, n_sexes, n_srv_fleets)) != 4)
      stop(paste("ObsSrvAgeComps is not the correct dimension. Should be n_regions, n_years, number of observed composition ages, n_sexes, n_srv_fleets"))
  }

  if(what == "ObsSrvLenComps") {
    if(sum(dim(x) == c(n_regions, n_years, n_lens, n_sexes, n_srv_fleets)) != 5)
      stop(paste("ObsSrvLenComps is not the correct dimension. Should be n_regions, n_years, n_lens, n_sexes, n_srv_fleets"))
  }

  if(what %in% c('ISS_SrvLenComps', 'ISS_SrvAgeComps')) {
    if(sum(dim(x) == c(n_regions, n_years, n_sexes, n_srv_fleets)) != 4)
      stop(paste(what, " is not the correct dimension. Should be n_regions, n_years, n_sexes, n_srv_fleets"))
  }

}

