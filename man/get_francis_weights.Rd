% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/do_francis.R
\name{get_francis_weights}
\alias{get_francis_weights}
\title{Computes Francis weights, which is used internally by do_francis_reweighting}
\usage{
get_francis_weights(
  n_regions,
  n_sexes,
  n_fleets,
  Use,
  ISS,
  Pred_array,
  Obs_array,
  weights,
  bins,
  comp_type
)
}
\arguments{
\item{n_regions}{Number of regions}

\item{n_sexes}{Number of sexes}

\item{n_fleets}{Number of fleets (fishery or survey)}

\item{Use}{Array from data list that specifies whether to use data that year}

\item{ISS}{Input sample size array}

\item{Pred_array}{Predicted values array dimensioned by n_regions, n_years, n_ages, n_sexes, n_fleets}

\item{Obs_array}{Observed values array dimensioned by n_regions, n_years, n_ages, n_sexes, n_fleets}

\item{weights}{Array of francis weights (NAs) to apply dimensioned by n_regions, n_years, n_sexes, n_fleets}

\item{bins}{Vector of bins used (age or length)}

\item{comp_type}{Matrix of composition structure types dimensioned by year and fleet}
}
\value{
Value for calculated francis weight
}
\description{
Computes Francis weights, which is used internally by do_francis_reweighting
}
\examples{
\dontrun{
# Function is used within do_francis_reweighting
get_francis_weights(n_regions, n_sexes, n_fleets, Use, ISS, Pred_arry, Obs_array, weights, bins, comp_type)
}
}
\keyword{internal}
