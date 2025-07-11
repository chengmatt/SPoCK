% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Get_Model_Fits.R
\name{get_idx_fits}
\alias{get_idx_fits}
\title{Gets index fits results}
\usage{
get_idx_fits(data, rep, year_labs)
}
\arguments{
\item{data}{Data list fed into RTMB}

\item{rep}{Report list output from RTMB}

\item{year_labs}{Year labels to use (vector)}
}
\value{
Fits to indices as a dataframe
}
\description{
Gets index fits results
}
\examples{
\dontrun{
idx_fits <- get_idx_fits(data = data, rep = rep, year_labs = seq(1960, 2024, 1))

idx_fits <- idx_fits \%>\%
  mutate(
    Idx = case_when(
      Type == "Fishery" & Year < 1995 ~ "Japanese Fishery CPUE Index",
      Type == "Fishery" & Year >= 1995 ~ "Domestic Fishery CPUE Index",
      Type == 'Survey' & Fleet == 1 ~ "Domestic LL Survey Relative Population Numbers",
      Type == 'Survey' & Fleet == 2 ~ "GOA Trawl Survey Biomass (kt)",
      Type == 'Survey' & Fleet == 3 ~ 'Japanese LL Survey Relative Population Numbers'
    )
  )
ggplot() +
  geom_line(idx_fits, mapping = aes(x = Year, y = value), lwd = 1.3, col = 'red') +
  geom_pointrange(idx_fits, mapping = aes(x = Year, y = obs, ymin = lci, ymax = uci), color = 'blue', pch = 1) +
  labs(x = "Year", y = 'Index') +
  theme_bw(base_size = 20) +
  facet_wrap(~Idx, scales = 'free', ncol = 2)

}
}
