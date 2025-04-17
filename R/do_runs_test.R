#' Runs test function taken from SS3 diags.
#'
#' @param x Vector of residuals
#' @param type Whether to use mean 0 assumption of mean of residuals (default = use mean 0)
#' @param mixing Type of test to do, less = left tailed test that detects positive autocorrelation, two.sided = two sided test that tests whether there is positive and/or negative autocorrealtion. The null is that there isn't any, rejecting the null (<0.05) indictes that there is some non-randomness.
#'
#' @returns List object with p value and limits for a three-sigma limit - (potential data outlier, where residual is > 3 standard deviations away from a mean of 0)
#' @export do_runs_test
#'
#' @import randtests
#' @examples
#' \dontrun{
#'  idx_fits <- get_idx_fits(data = data, rep = rep, year_labs = seq(1960, 2024, 1))
#'   idx_fits <- idx_fits %>%
#'     mutate(
#'       Idx = case_when(
#'         Type == "Fishery" & Year < 1995 ~ "Japanese Fishery CPUE Index",
#'         Type == "Fishery" & Year >= 1995 ~ "Domestic Fishery CPUE Index",
#'         Type == 'Survey' & Fleet == 1 ~ "Domestic LL Survey Relative Population Numbers",
#'         Type == 'Survey' & Fleet == 2 ~ "GOA Trawl Survey Biomass (kt)",
#'         Type == 'Survey' & Fleet == 3 ~ 'Japanese LL Survey Relative Population Numbers'
#'       )
#'     )
#'
#'   unique_idx <- unique(idx_fits$Idx)
#'   runs_all <- data.frame()
#'   for(i in 1:length(unique(idx_fits$Idx))) {
#'     tmp <- idx_fits %>% filter(Idx == unique_idx[i])
#'     runstest <- ssruns_sig3(x=as.numeric(tmp$resid),type="resid", mixing = "less")
#'     tmp_runs <- data.frame(p = runstest$p.runs, lwr = runstest$sig3lim[1], upr = runstest$sig3lim[2], Idx = unique_idx[i])
#'     runs_all <- rbind(runs_all, tmp_runs)
#'   } # end i
#'
#'   ggplot() +
#'     geom_point(idx_fits, mapping = aes(x = Year, y = resid)) +
#'     geom_segment(idx_fits, mapping = aes(x = Year, xend = Year, y = 0, yend = resid)) +
#'     geom_smooth(idx_fits, mapping = aes(x = Year, y = resid), se = F) +
#'     geom_hline(yintercept = 0, lty = 2) +
#'     geom_hline(runs_all, mapping = aes(yintercept = upr), lty = 2) +
#'     geom_hline(runs_all, mapping = aes(yintercept = lwr), lty = 2) +
#'     geom_text(data = runs_all, aes(x = -Inf, y = Inf, label = paste("p = ", round(p, 3))), hjust = -0.5, vjust = 8.2, size = 7)+
#'     labs(x = "Year", y = 'Residuals') +
#'     theme_bw(base_size = 20) +
#'     facet_wrap(~Idx, scales = 'free', ncol = 2)
#' }
do_runs_test <- function(x,
                         type=NULL,
                         mixing="two.sided"
                         ) {
  if(is.null(type)) type="resid"
  if(type=="resid"){
    mu = 0}else{mu = mean(x, na.rm = TRUE)}
  alternative=c("two.sided","left.sided")[which(c("two.sided", "less")%in%mixing)]
  # Average moving range
  mr  <- abs(diff(x - mu))
  amr <- mean(mr, na.rm = TRUE)
  # Upper limit for moving ranges
  ulmr <- 3.267 * amr
  # Remove moving ranges greater than ulmr and recalculate amr, Nelson 1982
  mr  <- mr[mr < ulmr]
  amr <- mean(mr, na.rm = TRUE)
  # Calculate standard deviation, Montgomery, 6.33
  stdev <- amr / 1.128
  # Calculate control limits
  lcl <- mu - 3 * stdev
  ucl <- mu + 3 * stdev
  if(nlevels(factor(sign(x)))>1){
    # Make the runs test non-parametric
    runstest = randtests::runs.test(x,threshold = 0,alternative = alternative)
    if(is.na(runstest$p.value)) p.value =0.001
    pvalue = round(runstest$p.value,3)} else {
      pvalue = 0.001
    }
  return(list(sig3lim=c(lcl,ucl),p.runs= pvalue))
}
