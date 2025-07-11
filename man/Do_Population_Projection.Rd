% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Do_Population_Projection.R
\name{Do_Population_Projection}
\alias{Do_Population_Projection}
\title{Do Population Projections}
\usage{
Do_Population_Projection(
  n_proj_yrs = 2,
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
  WAA_fish,
  MatAA,
  fish_sel,
  Movement,
  f_ref_pt = NULL,
  b_ref_pt = NULL,
  HCR_function = NULL,
  recruitment_opt = "inv_gauss",
  fmort_opt = "HCR",
  t_spawn,
  bh_rec_opt = NULL
)
}
\arguments{
\item{n_proj_yrs}{Number of projection years}

\item{n_regions}{Number of regions}

\item{n_ages}{Number of ages}

\item{n_sexes}{Number of sexes}

\item{sexratio}{Recruitment sex ratio}

\item{n_fish_fleets}{Number of fishery fleets}

\item{do_recruits_move}{Whether recruits move (0 == don't move, 1 == move)}

\item{recruitment}{Recruitment matrix dimensioned by n_regions, and n_yrs that we want to summarize across, or condition our projection on}

\item{terminal_NAA}{Terminal Numbers at Age dimensioned by n_regions, n_ages, n_sexes}

\item{terminal_F}{Terminal fishing mortality rate, dimensioned by n_regions, n_fish_fleets}

\item{natmort}{Natural mortality, dimensioned by n_regions, n_proj_yrs, n_ages, n_sexes}

\item{WAA}{Weight at age, dimensioned by n_regions, n_proj_yrs, n_ages, n_sexes}

\item{WAA_fish}{Weight at age for the fishery, dimensioned by n_regions, n_proj_yrs, n_ages, n_sexes, n_fish_fleets}

\item{MatAA}{Maturity at age, dimensioned by n_regions, n_proj_yrs, n_ages, n_sexes}

\item{fish_sel}{Fishery selectivity, dimensioned by n_regions, n_proj_yrs, n_ages, n_sexes, n_fish_fleets}

\item{Movement}{Movement, dimensioned by n_regions, n_regions, n_proj_yrs, n_ages, n_sexes}

\item{f_ref_pt}{Fishing mortality reference point dimensioned by n_regions and n_proj_yrs}

\item{b_ref_pt}{Biological reference point dimensioned by n_regions and n_proj_yrs}

\item{HCR_function}{Function describing a harvest control rule. The function should always have the following arguments: x, which represents SSB, frp, which takes inputs of fishery reference points, and brp, which takes inputs of biological reference points. Any additional arguments should be specified with defaults or hard coded / fixed within the function.}

\item{recruitment_opt}{Recruitment simulation option, where options are "inv_gauss", which simulates future recruitment based on the the recruitment values supplied using an inverse gaussian distribution, "mean_rec", which takes the mean of the recruitment values supplied for a given region, and "zero", which assumes that future recruitment does not occur}

\item{fmort_opt}{Fishing Mortality option, which includes "HCR", which modifies the F reference point using a user supplied HCR_function, or "Input", which uses projected F values supplied by the user.}

\item{t_spawn}{Fraction time of spawning used to compute projected SSB}

\item{bh_rec_opt}{A list object containing the following arguments:
\describe{
\item{recruitment_dd}{A value (0 or 1) indicating global (1) or local density dependence (0). In the case of a single region model, either local or global will give the same results}
\item{rec_lag}{A value indicating the number of years lagged that a given year's SSB produces recruits}
\item{R0}{The virgin recruitment parameter}
\item{Rec_Prop}{Recruitment apportionment values. In a single region model, this should be set at a value of 1. Dimensioned by n_regions}
\item{h}{Steepness values for the stock recruitment curve. Dimensioned by n_regions}
\item{WAA}{A weight-at-age array dimensioned by n_regions, n_ages, and n_sexes, where the reference year should utilize values from the first year}
\item{MatAA}{A maturity at age array dimensioned by n_regions, n_ages, and n_sexes, where the reference year should utilize values from the first year}
\item{natmort}{A natural mortality at age array dimensioned by n_regions, n_ages, and n_sexes, where the reference year should utilize values from the first year}
\item{SSB}{All SSB values estimated from a given model, dimensioned by n_regions and n_yrs}
}}
}
\value{
A list containing projected F, catch, SSB, and Numbers at Age. (Objects are generally dimensioned in the following order: n_regions, n_yrs, n_ages, n_sexes, n_fleets)
}
\description{
Do Population Projections
}
\examples{
\dontrun{
# Define HCR to use
HCR_function <- function(x, frp, brp, alpha = 0.05) {
  stock_status <- x / brp # define stock status
  # If stock status is > 1
  if(stock_status >= 1) f <- frp
  # If stock status is between brp and alpha
  if(stock_status > alpha && stock_status < 1) f <- frp * (stock_status - alpha) / (1 - alpha)
  # If stock status is less than alpha
  if(stock_status < alpha) f <- 0
  return(f)
}
rep <- obj$report(obj$env$last.par.best) # need to have an RTMB object first
# Setup necessary inputs
n_sims <- 1000
t_spawn <- 0
sexratio <- 0.5
n_proj_yrs <- 15
n_regions <- 1
n_ages <- length(data$ages)
n_sexes <- data$n_sexes
n_fish_fleets <- 2
do_recruits_move <- 0
terminal_NAA <- array(obj$rep$NAA[,length(data$years),,], dim = c(n_regions, n_ages, n_sexes))
WAA <- array(rep(data$WAA[,length(data$years),,], each = n_proj_yrs), dim = c(n_regions, n_proj_yrs, n_ages, n_sexes)) # weight at age
WAA_fish <- array(rep(data$WAA_fish[,length(data$years),,,], each = n_proj_yrs), dim = c(n_regions, n_proj_yrs, n_ages, n_sexes, n_fish_fleets)) # weight at age for fishery
MatAA <- array(rep(data$MatAA[,length(data$years),,], each = n_proj_yrs), dim = c(n_regions, n_proj_yrs, n_ages, n_sexes)) # maturity at age
fish_sel <- array(rep(obj$rep$fish_sel[,length(data$years),,,], each = n_proj_yrs), dim = c(n_regions, n_proj_yrs, n_ages, n_sexes, n_fish_fleets)) # selectivity
Movement <- array(rep(obj$rep$Movement[,,length(data$years),,], each = n_proj_yrs), dim = c(n_regions, n_regions, n_proj_yrs, n_ages, n_sexes))
terminal_F <- array(obj$rep$Fmort[,length(data$years),], dim = c(n_regions, n_fish_fleets))
natmort <- array(obj$rep$natmort[,length(data$years),,], dim = c(n_regions, n_proj_yrs, n_ages, n_sexes))
recruitment <- array(obj$rep$Rec[,20:(length(data$years) - 2)], dim = c(n_regions, length(20:length(data$years) - 2)))

# Define reference points
spr_35 <- Get_Reference_Points(data = data,
                               rep = rep,
                               SPR_x = 0.35, t_spwn = 0, sex_ratio_f = 0.5,
                               calc_rec_st_yr = 20, rec_age = 2)

spr_40 <- Get_Reference_Points(data = data,
                               rep = rep,
                               SPR_x = 0.4, t_spwn = 0, sex_ratio_f = 0.5,
                               calc_rec_st_yr = 20, rec_age = 2)

spr_60 <- Get_Reference_Points(data = data,
                               rep = rep,
                               SPR_x = 0.6, t_spwn = 0, sex_ratio_f = 0.5,
                               calc_rec_st_yr = 20, rec_age = 2)

# Extract reference points
b40 <- spr_40$b_ref_pt
b60 <- spr_60$b_ref_pt
f40 <- spr_40$f_ref_pt
f35 <- spr_35$f_ref_pt
f60 <- spr_60$f_ref_pt
# Define the F used for each scenario (Based on BSAI Intro Report)
proj_inputs <- list(
  # Scenario 1 - Using HCR to adjust maxFABC
  list(f_ref_pt = array(f40, dim = c(n_regions, n_proj_yrs)),
       b_ref_pt = array(b40, dim = c(n_regions, n_proj_yrs)),
       fmort_opt = 'HCR'
  ),
  # Scenario 2 - Using HCR to adjust maxFABC based on last year's value (constant fraction - author specified F)
  list(f_ref_pt = array(f40 * (f40 / 0.086), dim = c(n_regions, n_proj_yrs)),
       b_ref_pt = array(b40, dim = c(n_regions, n_proj_yrs)),
       fmort_opt = 'HCR'
  ),
  # Scenario 3 - Using an F input of last 5 years average F, and
  list(f_ref_pt = array(mean(rowSums(sabie_rtmb_model$rep$Fmort[1, 60:64, ])), dim = c(n_regions, n_proj_yrs)),
       b_ref_pt = NULL,
       fmort_opt = 'Input'
  ),
  # Scenario 4 - Using HCR to adjust F60
  list(f_ref_pt = array(f60, dim = c(n_regions, n_proj_yrs)),
       b_ref_pt = array(b40, dim = c(n_regions, n_proj_yrs)),
       fmort_opt = 'HCR'
  ),
  # Scenario 5 - F is set at 0
  list(f_ref_pt = array(0, dim = c(n_regions, n_proj_yrs)),
       b_ref_pt = NULL,
       fmort_opt = 'Input'
  ),
  # Scenario 6 - Using HCR to adjust FOFL
  list(f_ref_pt = array(f35, dim = c(n_regions, n_proj_yrs)),
       b_ref_pt = array(b40, dim = c(n_regions, n_proj_yrs)),
       fmort_opt = 'HCR'
  ),
  # Scenario 7 - Using HCR to adjust FABC in first 2 projection years, and then later years are adjusting FOFL
  list(f_ref_pt = array(c(rep(f40, 2), rep(f35, n_proj_yrs - 2)), dim = c(n_regions, n_proj_yrs)),
       b_ref_pt = array(b40, dim = c(n_regions, n_proj_yrs)),
       fmort_opt = 'HCR'
  )
)

# store outputs
all_scenarios_f <- array(0, dim = c(n_regions, n_proj_yrs, n_sims, length(proj_inputs)))
all_scenarios_ssb <- array(0, dim = c(n_regions, n_proj_yrs, n_sims, length(proj_inputs)))
all_scenarios_catch <- array(0, dim = c(n_regions, n_proj_yrs, n_fish_fleets, n_sims, length(proj_inputs)))

for (i in seq_along(proj_inputs)) {
  for (sim in 1:n_sims) {

    # do population projection
    out <- Do_Population_Projection(n_proj_yrs = n_proj_yrs,
                                    n_regions = n_regions,
                                    n_ages = n_ages,
                                    n_sexes = n_sexes,
                                    sexratio = sexratio,
                                    n_fish_fleets = n_fish_fleets,
                                    do_recruits_move = do_recruits_move,
                                    recruitment = recruitment,
                                    terminal_NAA = terminal_NAA,
                                    terminal_F = terminal_F,
                                    natmort = natmort,
                                    WAA = WAA,
                                    WAA_fish = WAA_fish,
                                    MatAA = MatAA,
                                    fish_sel = fish_sel,
                                    Movement = Movement,
                                    f_ref_pt = proj_inputs[[i]]$f_ref_pt,
                                    b_ref_pt = proj_inputs[[i]]$b_ref_pt,
                                    HCR_function = HCR_function,
                                    recruitment_opt = "inv_gauss",
                                    fmort_opt = proj_inputs[[i]]$fmort_opt,
                                    t_spawn = t_spawn
    )

    all_scenarios_ssb[,,sim,i] <- out$proj_SSB
    all_scenarios_catch[,,,sim,i] <- out$proj_Catch
    all_scenarios_f[,,sim,i] <- out$proj_F[,-(n_proj_yrs+1)] # remove last year, since it's not used
  } # end sim loop
  print(i)
} # end i loop

# If users were to specify "bh_rec" for recruitment_opt, a list of specifications for projecting deterministic recruitment is required. An example
# of this is provided below:
bh_rec_opt <- list(
  recruitment_dd = 1,
  rec_lag = 1,
  R0 = rep$R0,
  h = rep$h_trans,
  Rec_Prop = 1,
  WAA = array(data$WAA[,1,,], dim = c(1, n_ages, n_sexes)),
  MatAA = array(data$MatAA[,1,,], dim = c(1, n_ages, n_sexes)),
  natmort = array(data$Fixed_natmort[,1,,], dim = c(1, n_ages, n_sexes)),
  SSB = rep$SSB
)
}
}
