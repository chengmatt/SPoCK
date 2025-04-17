#' Get Selectivity
#'
#' @param Selex_Model integer for selectivity model
#' @param ln_Pars log selectivity parameters
#' @param Age vector of Ages
#' @param TimeVary_Model Time-varying selectivity model == 0, constant or blocked, == 1 iid, == 2 random walk, == 3 3dar1 marginal, == 4 3dar1 conditional, ==5 2dar1
#' @param ln_seldevs selectivity deviations (either in array or matrix and needs to match the time_vary model specification)
#' @param Year Year index value
#' @param Region Region index value
#' @param Sex Sex index value
#'
#' @export Get_Selex
Get_Selex = function(Selex_Model, TimeVary_Model, ln_Pars, ln_seldevs, Region, Year, Age, Sex) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  selex = rep(0, length(Age)) # Temporary container vector
  if(Selex_Model == 0) { # logistic selectivity (a50 and slope)
    # Extract out and exponentiate the parameters here
    a50 = exp(ln_Pars[1]); # a50
    k = exp(ln_Pars[2]); # slope
    if(TimeVary_Model %in% c(1,2)) {
      a50 = a50 * exp(ln_seldevs[Region, Year, 1, Sex, 1]) # a50 parameter varying
      k = k * exp(ln_seldevs[Region, Year, 2, Sex, 1]) # slope parameter varying
    } # end if iid or random walk

    selex = 1 / (1 + exp(-k * (Age - a50))) # return parmetric form

    # 3dgmrf model or 2dar1 (sel devs dimensioned as region, year, age, sex
    if(TimeVary_Model %in% c(3,4,5)) selex = selex * exp(ln_seldevs[Region,Year,,Sex, 1]) # varies semi-parametriclly
  }

  if(Selex_Model == 1) { # gamma dome-shaped selectivity
    # Extract out and exponentiate the parameters here
    amax = exp(ln_Pars[1]) # age at max selex
    delta = exp(ln_Pars[2]) # slope parameter
    if(TimeVary_Model %in% c(1,2)) {
      amax = amax * exp(ln_seldevs[Region, Year, 1, Sex, 1]) # amax parameter varying
      delta = delta * exp(ln_seldevs[Region, Year, 2, Sex, 1]) # delta parameter varying
    } # end if iid or random walk
    # Now, calculate/derive power parameter + selex values
    p = 0.5 * (sqrt( amax^2 + (4 * delta^2)) - amax)

    selex = (Age / amax)^(amax/p) * exp( (amax - Age) / p ) # return parametric form

    # 3dgmrf model or 2dar1 (sel devs dimensioned as region, year, age, sex
    if(TimeVary_Model %in% c(3,4,5)) selex = selex * exp(ln_seldevs[Region,Year,,Sex, 1]) # varies semi-parametriclly
  }

  if(Selex_Model == 2) { # power function selectivity
    # Extract out and exponentiate the parameters here
    power = exp(ln_Pars[1]); # power parameter
    if(TimeVary_Model %in% c(1,2)) {
      power = power * exp(ln_seldevs[Region, Year, 1, Sex, 1]) # power parameter varying
    } # end if iid or random walk

    selex = 1 / Age^power # return parametric form

    # 3dgmrf model or 2dar1 (sel devs dimensioned as region, year, age, sex
    if(TimeVary_Model %in% c(3,4,5)) selex = selex * exp(ln_seldevs[Region,Year,,Sex, 1]) # varies semi-parametriclly
  }

  if(Selex_Model == 3) { # logistic selectivity (a50 and a95)
    # Extract out and exponentiate the parameters here
    a50 = exp(ln_Pars[1]); # a50
    a95 = exp(ln_Pars[2]); # a95

    if(TimeVary_Model %in% c(1,2)) {
      a50 = a50 * exp(ln_seldevs[Region, Year, 1, Sex, 1]) # a50 parameter varying
      a95 = a95 * exp(ln_seldevs[Region, Year, 2, Sex, 1]) # a95 parameter varying
    } # end if iid or random walk

    selex = 1 / (1+19^((a50-Age)/a95)) # 19 b/c 0.95 / (1 - 0.95) return parametric form

    # 3dgmrf model or 2dar1 (sel devs dimensioned as region, year, age, sex
    if(TimeVary_Model %in% c(3,4,5)) selex = selex * exp(ln_seldevs[Region,Year,,Sex, 1]) # varies semi-parametriclly
  }

  return(selex)
} # end function

