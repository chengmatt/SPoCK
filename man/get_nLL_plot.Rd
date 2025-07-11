% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Make_Plots_Tables.R
\name{get_nLL_plot}
\alias{get_nLL_plot}
\title{Get plot of negative log likelihood values}
\usage{
get_nLL_plot(rep, model_names)
}
\arguments{
\item{rep}{List of n_models of `SPoRC` report lists}

\item{model_names}{Vector of model names}
}
\value{
Plot and tables of negative log likelihood values across models
}
\description{
Get plot of negative log likelihood values
}
\examples{
\dontrun{
get_nLL_plot(list(rep1, rep2), c("Model1", "Model2"))
}
}
