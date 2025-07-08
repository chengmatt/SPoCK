#' Calculate Selectivity
#'
#' Computes selectivity using one of several parametric or semi-parametric models.
#' Supports both constant and time-varying selectivity, including random effects and GMRF-based deviations.
#'
#' @param Selex_Model Integer specifying the selectivity model:
#'   \describe{
#'     \item{0}{Logistic selectivity: uses b50 and slope parameters}
#'     \item{1}{Gamma-shaped (dome) selectivity: uses bin-at-peak and delta parameters}
#'     \item{2}{Power function selectivity: decreasing selectivity with bin}
#'     \item{3}{Logistic selectivity using b50 and b95}
#'     \item{4}{Double-normal (dome-shaped) selectivity with plateau and flexible tails}
#'   }
#'
#' @param TimeVary_Model Integer specifying time variation structure:
#'   \describe{
#'     \item{0}{No time variation (constant or blocked)}
#'     \item{1}{IID deviations}
#'     \item{2}{Random walk over time}
#'     \item{3}{3D AR1-GMRF marginal}
#'     \item{4}{3D AR1-GMRF conditional}
#'     \item{5}{2D AR1-GMRF}
#'   }
#'
#' @param ln_Pars Vector of log-transformed selectivity parameters. Interpretation depends on `Selex_Model`.
#' @param ln_seldevs Array of selectivity deviations (may be log-scale), dimensioned as:
#'   [n_regions, n_years, n_bins, n_sexes, 1]. Used for time-varying or semi-parametric selectivity.
#' @param Region Integer index for region
#' @param Year Integer index for year
#' @param Bin Numeric vector of bins to compute selectivity for
#' @param Sex Integer index for sex
#'
#' @return A numeric vector of selectivity values corresponding to the bins specified in the model.
#'
#' @details
#' Selectivity parameters are transformed internally (typically using `exp()` or logistic transformations)
#' to ensure they remain in valid ranges. Deviations (`ln_seldevs`) apply multiplicatively to these transformed
#' parameters when time-varying models are used. For semi-parametric models (TimeVary_Model 3–5), deviations are applied directly to the resulting selectivity curve.
#'
#' @keywords internal
Get_Selex = function(Selex_Model,
                     TimeVary_Model,
                     ln_Pars,
                     ln_seldevs,
                     Region,
                     Year,
                     Bin,
                     Sex) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  selex = rep(0, length(Bin)) # Temporary container vector

  if(Selex_Model == 0) { # logistic selectivity (b50 and slope)
    # Extract out and exponentiate the parameters here
    b50 = exp(ln_Pars[1]); # b50
    k = exp(ln_Pars[2]); # slope

    if(TimeVary_Model %in% c(1:2)) {
      b50 = b50 * exp(ln_seldevs[Region, Year, 1, Sex, 1]) # b50 parameter varying
      k = k * exp(ln_seldevs[Region, Year, 2, Sex, 1]) # slope parameter varying
    } # end if iid or random walk

    selex = 1 / (1 + exp(-k * (Bin - b50))) # return parmetric form
  }

  if(Selex_Model == 1) { # gamma dome-shaped selectivity
    # Extract out and exponentiate the parameters here
    bmax = exp(ln_Pars[1]) # Bin at max selex
    delta = exp(ln_Pars[2]) # slope parameter

    if(TimeVary_Model %in% c(1:2)) {
      bmax = bmax * exp(ln_seldevs[Region, Year, 1, Sex, 1]) # bmax parameter varying
      delta = delta * exp(ln_seldevs[Region, Year, 2, Sex, 1]) # delta parameter varying
    } # end if iid or random walk

    # Now, calculate/derive power parameter + selex values
    p = 0.5 * (sqrt( bmax^2 + (4 * delta^2)) - bmax)
    selex = (Bin / bmax)^(bmax/p) * exp( (bmax - Bin) / p ) # return parametric form
  }

  if(Selex_Model == 2) { # power function selectivity
    # Extract out and exponentiate the parameters here
    power = exp(ln_Pars[1]); # power parameter

    if(TimeVary_Model %in% c(1:2)) {
      power = power * exp(ln_seldevs[Region, Year, 1, Sex, 1]) # power parameter varying
    } # end if iid or random walk

    selex = 1 / Bin^power # return parametric form
  }

  if(Selex_Model == 3) { # logistic selectivity (b50 and b95)

    # Extract out and exponentiate the parameters here
    b50 = exp(ln_Pars[1]); # b50
    b95 = exp(ln_Pars[2]); # b95

    if(TimeVary_Model %in% c(1:2)) {
      b50 = b50 * exp(ln_seldevs[Region, Year, 1, Sex, 1]) # b50 parameter varying
      b95 = b95 * exp(ln_seldevs[Region, Year, 2, Sex, 1]) # b95 parameter varying
    } # end if iid or random walk

    selex = 1 / (1+19^((b50-Bin)/b95)) # 19 b/c 0.95 / (1 - 0.95) return parametric form
  }

  if(Selex_Model == 4) {

    # define bin ranges for double normal here
    max_x_val <- max(Bin)
    midbin <- Bin

    # Extract and transform parameters here
    p1trans <- min(Bin) + (max(Bin) - min(Bin)) * RTMB::plogis(ln_Pars[1]) # peak bin at plateau
    p2trans <- p1trans + 1 + (0.99 + max_x_val - p1trans - 1)/(1 + exp(-1.0 * ln_Pars[2])) # width of plateau
    p3trans <- exp(ln_Pars[3]) # ascending width
    p4trans <- exp(ln_Pars[4]) # descending width
    p5trans <- 1/(1 + exp(-1.0 * ln_Pars[5])) # selectivity at first bin
    p6trans <- 1/(1 + exp(-1.0 * ln_Pars[6])) # selectivity at last bin

    if(TimeVary_Model %in% c(1:2)) {
      p1trans = p1trans * exp(ln_seldevs[Region, Year, 1, Sex, 1]) # p1 parameter varying
      p2trans = p2trans * exp(ln_seldevs[Region, Year, 2, Sex, 1]) # p2 parameter varying
      p3trans = p3trans * exp(ln_seldevs[Region, Year, 3, Sex, 1]) # p3 parameter varying
      p4trans = p4trans * exp(ln_seldevs[Region, Year, 4, Sex, 1]) # p4 parameter varying
      p5trans = p5trans * exp(ln_seldevs[Region, Year, 5, Sex, 1]) # p5 parameter varying
      p6trans = p6trans * exp(ln_seldevs[Region, Year, 6, Sex, 1]) # p6 parameter varying
    } # end if iid or random walk

    # construct selectivity function
    asc <- exp(-((midbin - p1trans)^2/p3trans))
    asc.scaled <- (p5trans + (1 - p5trans) * (asc - 0)/(1 - 0))
    desc <- exp(-((midbin - p2trans)^2/p4trans))
    stj <- exp(-((40 - p2trans)^2/p4trans))
    des.scaled <- (1 + (p6trans - 1) * (desc - 1) /(stj - 1))
    join1 <- 1/(1 + exp(-(20 * (midbin - p1trans)/(1 + abs(midbin - p1trans))))) # joiner functions
    join2 <- 1/(1 + exp(-(20 * (midbin - p2trans)/(1 + abs(midbin - p2trans))))) # joiner functions
    selex <- asc.scaled * (1 - join1) + join1 * (1 * (1 - join2) + des.scaled * join2) # return parameteric form
    selex[1] <- p5trans # return parameteric form
  }

  # 3dgmrf model or 2dar1 (sel devs dimensioned as region, year, bin, sex)
  if(TimeVary_Model %in% c(3:5)) selex = selex * exp(ln_seldevs[Region,Year,,Sex, 1]) # varies semi-parametriclly

  return(selex)
} # end function

