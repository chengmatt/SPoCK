% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utility_Functions.R
\name{bisection_F}
\alias{bisection_F}
\title{Go from TAC to Fishing Mortality using bisection}
\usage{
bisection_F(
  f_guess,
  catch,
  NAA,
  WAA,
  natmort,
  fish_sel,
  n.iter = 20,
  lb = 0,
  ub = 2
)
}
\arguments{
\item{f_guess}{Initial guess of F}

\item{catch}{Provided catch values}

\item{NAA}{Numbers, dimensioned by ages, and sexes}

\item{WAA}{Weight, dimensioned by ages and sexes}

\item{natmort}{Natural mortality dimensioned by ages and sex}

\item{fish_sel}{Fishery selectivity, dimesnioned by ages and sex}

\item{n.iter}{Number of iterations for bisection}

\item{lb}{Lower bound of F}

\item{ub}{Upper bound of F}
}
\value{
Fishing mortality values
}
\description{
Go from TAC to Fishing Mortality using bisection
}
