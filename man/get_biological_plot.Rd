% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Make_Plots_Tables.R
\name{get_biological_plot}
\alias{get_biological_plot}
\title{Get Plots of Biological Quantities}
\usage{
get_biological_plot(data, rep, model_names)
}
\arguments{
\item{data}{List of n_models of `SPoRC` data lists}

\item{rep}{List of n_models of `SPoRC` report lists}

\item{model_names}{Vector of model names}
}
\value{
A list of plots for terminal year movement, natural mortality, weight-at-age, and maturity at age across models
}
\description{
Get Plots of Biological Quantities
}
\examples{
\dontrun{
get_biological_plot(list(data1, data2), list(rep1, rep2), c("Model1", "Model2"))
}
}
