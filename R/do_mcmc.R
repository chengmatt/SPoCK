#' Run MCMC using rtmbstan
#'
#' @param obj Object built from ADFUN via RTMB
#' @param chains Number of chains to run
#' @param cores Number of cores to use
#' @param iter Number of MCMC iterations to run
#' @param thin Thinning rate
#' @param bounds Uniform bounds to constrain parameter bounds
#' @param laplace Whether or not to do laplace approximation for random effects and MCMC for fixed effects
#' @param ... Additional arguments for tmbstan (e.g., adapt_delta)
#'
#' @import tmbstan
#' @import RTMB
#' @importFrom stats rnorm nlminb optimHess
#' @returns MCMC list object from rtmbstan
#' @export do_mcmc
#'
#' @examples
#' \dontrun{
#' obj <- RTMB::MakeADFun(cmb(SPoRC_rtmb, data), parameters = parameters, map = mapping)
#' mcmc <- do_mcmc(obj = obj,
#'                 chains = 4,
#'                 cores = 4,
#'                 iter = 10000,
#'                 thin = 10,
#'                 bounds = 1.5,
#'                 laplace = FALSE,
#'                 adapt_delta = 0.99)
#' saveRDS(mcmc, here('output', 'MCMC_Model_23.5.rds'))
#' }
do_mcmc <- function(obj,
                    chains,
                    cores,
                    iter,
                    thin,
                    bounds = Inf,
                    laplace = FALSE,
                    ...
                    ) {

  # Optimize MLE for use in tuning MCMC
  optim <- stats::nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 1e6, eval.max = 1e6, rel.tol = 1e-15))

  # run newton steps to get lower
  tryCatch(expr =
             for(i in 1:3) {
               g = as.numeric(obj$gr(optim$par))
               h = optimHess(optim$par, fn = obj$fn, gr = obj$gr)
               optim$par = optim$par - solve(h,g)
               optim$objective = obj$fn(optim$par)
             }
           , error = function(e){e}, warning = function(w){w})

  # Setup MCMC cores
  options(mc.cores = cores)

  # Define bounds
  upper <- optim$par + bounds
  lower <- optim$par - bounds

  # Run MCMC
  fit <- tmbstan::tmbstan(obj = obj, iter = iter, thin = thin, chains = chains,
                          upper = upper, laplace = laplace,
                          lower = lower, control = list(...),
                          init = optim$par, open_progress = FALSE)

  return(fit)
} # end function

#' Gets parameter names from mcmc object
#'
#' @param mcmc_obj MCMC object from rtmbstan
#'
#' @importFrom rstan extract
#' @returns Returns a character vector of parameter names
#' @export mcmc_par_names
#'
#' @examples
#' \dontrun{
#' par_names <- mcmc_par_names(mcmc_obj)
#' }
mcmc_par_names <- function(mcmc_obj) {
  vals <- as.data.frame(rstan::extract(mcmc_obj))
  return(colnames(vals))
}


#' Model diagnostics for MCMC
#'
#' @param mcmc_obj MCMC object from tmbstan
#' @param tmb_sdrep SD report fomr TMB
#' @param plot Whether or not to plot outputs
#' @param pars Parameters you want to plot out
#' @importFrom rstan extract check_hmc_diagnostics summary
#' @importFrom GGally ggpairs
#' @importFrom bayesplot mcmc_trace
#' @importFrom cowplot plot_grid
#' @importFrom tidyr drop_na
#'
#' @returns A list of plots with traditional MCMC diagnostics
#' @export mcmc_diag
#'
#' @examples
#' \dontrun{
#' diagnostic_plots <- mcmc_diag(mcmc_obj, tmb_sdrep, TRUE, pars = 'ln_M')
#' }
mcmc_diag <- function(mcmc_obj,
                      tmb_sdrep,
                      plot = TRUE,
                      pars = NULL
                      ) {

  # MLE summary
  mle_sd <- data.frame(summary(tmb_sdrep)) %>% dplyr::mutate(lwr_95 = Estimate - 1.96 * Std..Error, upr_95 = Estimate + 1.96 * Std..Error) %>% dplyr::rename(mle_se = Std..Error, mle_est = Estimate)
  mle_sd$par <- ave(as.character(rownames(summary(tmb_sdrep))) , as.character(rownames(summary(tmb_sdrep))), FUN = function(z) if (length(z) == 1) z else paste0(z, ".", seq_along(z))) # make unique names

  # MCMC summary
  mcmc_summary <- rstan::summary(mcmc_obj)
  vals <- as.data.frame(rstan::extract(mcmc_obj))
  vals_long <- vals %>% pivot_longer(names_to = 'par', values_to = 'val', cols = everything()) # convert to long format

  # Get Rhats and ess
  diag_df <- data.frame(parnames = rownames(mcmc_summary$summary),
                        n_eff = mcmc_summary$summary[,'n_eff'],
                        R_hat = mcmc_summary$summary[,'Rhat'])

  if(any(diag_df$n_eff < 400)){
    message(paste("Some parameters have an effective sample size < 400. Inference may not be reliable. Consider reparameterizing or thinning more samples.",
                  "These parameters are:", diag_df$parnames[which(diag_df$n_eff < 400)]))
  } else message("All parameters have an effective sample size > 400. ")

  if(!is.na(any(diag_df$R_hat > 1.01))) {
    if(any(diag_df$R_hat > 1.01)) {
      message(paste("Some parameters have an Rhat > 1.01, which suggests that MCMC chains are not well-mixed and inference may not be reliable. Consider reparameterizing, increasing adapt_delta, or increasing the number of iterations.",
                    "These parameters are: ", diag_df$parnames[which(diag_df$R_hat > 1.01)]))
    } else message("All parameters have an Rhat < 1.01")
  } else {
    message("Some Rhats are NAs. Model is not converged!")
  }

  rstan::check_hmc_diagnostics(mcmc_obj)

  if(plot == TRUE) {

    # if no parameters specified
    if(is.null(pars)) {
      message("No parameters specified. Plotting the 10 least well-mixed parameters.")

      # find the 10 parameters that are least well mixed
      least_mixed_10 <- order(diag_df$R_hat, decreasing = TRUE)[1:10]
      least_mixed_10_df <- vals[,least_mixed_10]
      pars <- colnames(least_mixed_10_df)

      # reformat parameter names for trace
      base_name <- sub("\\.\\d+$", "", pars)
      number <- as.numeric(sub(".*\\.", "", pars))
      par_trace <- c(paste0(base_name[!is.na(number)], "[", number[!is.na(number)], "]"), paste(base_name[is.na(number)]))

      # Pairs plot
      pairs_plot <- GGally::ggpairs(least_mixed_10_df) + theme_bw(base_size = 15)

      # Trace plot
      trace_plot <- bayesplot::mcmc_trace(mcmc_obj, pars = par_trace) + theme_bw(base_size = 15)

      # Density plot with overlapping MLE
      dens_plot <- ggplot() +
        geom_density(vals_long %>% filter(par %in% pars), mapping = aes(x = val), fill = 'grey', color = 'black') +
        geom_pointrange(mle_sd %>% filter(par %in% pars), mapping =  aes(x = mle_est, y = 0, xmin = lwr_95, xmax = upr_95), fatten = 8, linewidth = 1.3) +
        facet_wrap(~par, scales = 'free') +
        labs(x = 'Parmeter Value', y = 'Density') +
        theme_bw(base_size = 15)
    }

    # parameters specified
    if(!is.null(pars)) {
      pars_df <- vals[,pars] # extract parameters we want to see

      # reformat parameter names for trace
      base_name <- sub("\\.\\d+$", "", pars)
      number <- as.numeric(sub(".*\\.", "", pars))
      par_trace <- c(paste0(base_name[!is.na(number)], "[", number[!is.na(number)], "]"), paste(base_name[is.na(number)]))
      par_trace <- par_trace[par_trace != '[]']

      # Pairs plot
      if(length(pars) == 1) pairs_plot <- NULL
      else pairs_plot <- ggpairs(pars_df) + theme_bw(base_size = 15)

      # Trace plot
      trace_plot <- bayesplot::mcmc_trace(mcmc_obj, pars = par_trace) + theme_bw(base_size = 15)

      # Density plot with overlapping MLE
      dens_plot <- ggplot() +
        geom_density(vals_long %>% filter(par %in% pars), mapping = aes(x = val), fill = 'grey', color = 'black') +
        geom_pointrange(mle_sd %>% filter(par %in% pars), mapping =  aes(x = mle_est, y = 0, xmin = lwr_95, xmax = upr_95), fatten = 8, linewidth = 1.3) +
        facet_wrap(~par, scales = 'free') +
        labs(x = 'Parmeter Value', y = 'Density') +
        theme_bw(base_size = 15)
    }

    # Compare mcmc and mle
    mcmc_mle_df <- vals_long %>% dplyr::group_by(par) %>% dplyr::summarize(mcmc_mean = median(val), mcmc_sd = sd(val)) %>%
      dplyr::left_join(mle_sd %>% select(par, mle_est, mle_se), by = 'par') %>% tidyr::drop_na()

    mean_plot <- ggplot(mcmc_mle_df, aes(x = mcmc_mean, y = mle_est)) +
      geom_point(size = 3, alpha = 0.5) +
      geom_abline(slope = 1, intercept = 0, lty = 2) +
      labs(x = "MCMC Median Estimate", y = "MLE Estimate") +
      theme_bw(base_size = 15)

    se_plot <- ggplot(mcmc_mle_df, aes(x = mcmc_sd, y = mle_se)) +
      geom_point(size = 3, alpha = 0.5) +
      geom_abline(slope = 1, intercept = 0, lty = 2) +
      labs(x = "MCMC SE", y = "MLE SE") +
      theme_bw(base_size = 15)

    comparision_mle_mcmc <- cowplot::plot_grid(mean_plot, se_plot)

  }

  return(list(pairs_plot, trace_plot, dens_plot, comparision_mle_mcmc))
}


#' Derive SSB and Recruitment from MCMC using posterior samples
#'
#' @param mcmc_obj MCMC object from tmbstan
#' @param tmb_obj RTMB object from MakeADFUN
#'
#' @import doSNOW
#' @import parallel
#' @importFrom data.table rbindlist
#' @returns returns a dataframe of recruitment and ssb posterior samples and a plot of ssb and recruitment with mean, lwr 95 quantile, and upr 95 quantile.
#' @export get_mcmc_ssb_rec
#'
get_mcmc_ssb_rec <- function(mcmc_obj,
                             tmb_obj
                             ) {

  # set up cores to extract MCMC results
  ncores <- detectCores()
  cl <- makeCluster(ncores - 2)
  registerDoSNOW(cl)

  # extract values
  posterior <- as.matrix(mcmc_obj) # posteriors
  n_iters <- nrow(posterior) # number of iterations to loop through

  # set up progress bar
  pb <- txtProgressBar(max = n_iters, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  # Get posteriors for derived variables
  derived_post <- foreach(iter = 1:n_iters, .packages = c("RTMB"), .options.snow = opts) %dopar% {

    # set up parameter data frame for posterior to join to so the parameter ordering is the same
    par_dat <- data.frame(tmb_obj$par, names = names(tmb_obj$par))
    par_dat$names <- ave(par_dat$names, par_dat$names, FUN = function(z) {
      if (length(z) == 1) return(z)
      paste0(z[1], "[", seq_along(z), "]")
    }) # making parameter names to be consistent with MCMC posterior naming format
    tmp_post_df <- data.frame(mcmc_pars = posterior[iter,-ncol(posterior)], names = colnames(posterior)[-ncol(posterior)])
    par_dat <- par_dat %>% dplyr::left_join(tmp_post_df, by = 'names') # join mcmc samples to correct parameter ordering

    post_rep <- tmb_obj$report(par_dat$mcmc_pars) # get derived variables here

    ssb_rec_post_df <- reshape2::melt(post_rep$SSB) %>% dplyr::rename(Region = Var1, Year = Var2) %>% dplyr::mutate(type = 'SSB') %>%
      dplyr::bind_rows(reshape2::melt(post_rep$Rec) %>% dplyr::rename(Region = Var1, Year = Var2) %>% dplyr::mutate(type = 'Rec')) %>%
      dplyr::mutate(iter = iter)

  } # end for each loop

  derived_post_df <- data.table::rbindlist(derived_post) # get into df

  # Summarize recruitment and ssb quantiles
  summarized <- derived_post_df %>%
    dplyr::mutate(Year = Year + 1959) %>%  # Adding years back
    dplyr::group_by(Year, type) %>%
    dplyr::summarize(mean_value = mean(value),
              lwr_95 = quantile(value, 0.025),
              upr_95 = quantile(value, 0.975))

  # plot
  summary_plot <- ggplot(summarized, aes(x = Year, y = mean_value, ymin = lwr_95, ymax = upr_95)) +
    geom_line(lty = 2) +
    geom_ribbon(alpha = 0.3) +
    facet_wrap(~type, scales = 'free') +
    theme_bw(base_size = 20) +
    labs(x = 'Year', y = 'Value')

  return(list(derived_post_df, summary_plot))

}
