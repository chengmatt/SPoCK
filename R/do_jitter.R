#' Function to conduct jitter analysis
#'
#' @param data Data list to make obj
#' @param parameters Parameter list to make obj
#' @param mapping Mapping list to make obj
#' @param random Character of random effects
#' @param sd sd for jitter (additive)
#' @param n_jitter Number of jitters to do
#' @param n_newton_loops Number of newton loops to do
#' @param do_par Whether to do paralleizaiton or not (boolean)
#' @param n_cores Number of cores to use
#'
#' @import RTMB
#' @import future.apply
#' @import progressr
#' @importFrom reshape2 melt
#' @import dplyr
#' @importFrom stats rnorm nlminb
#' @returns Dataframe of jitter values
#' @export do_jitter
#'
#' @examples
#' \dontrun{
#'    library(ggplot2)
#'    # get jitter values
#'    jit <- do_jitter(data = data,
#'                  parameters = parameters,
#'                  mapping = mapping,
#'                  random = NULL,
#'                  sd = 0.1,
#'                  n_jitter = 100,
#'                  n_newton_loops = 3,
#'                  do_par = TRUE,
#'                  n_cores = 8)
#'
#'    # get proportion converged
#'    prop_converged <- jit %>%
#'    filter(Year == 1, Type == 'Recruitment') %>%
#'      summarize(prop_conv = sum(Hessian) / length(Hessian))
#'
#'    # get final model results
#'    final_mod <- reshape2::melt(sabie_rtmb_model$rep$SSB) %>% rename(Region = Var1, Year = Var2) %>%
#'    mutate(Type = 'SSB') %>%
#'      bind_rows(reshape2::melt(sabie_rtmb_model$rep$Rec) %>%
#'      rename(Region = Var1, Year = Var2) %>% mutate(Type = 'Recruitment'))
#'
#'    # comparison of SSB and recruitment
#'   ggplot() +
#'     geom_line(jit, mapping = aes(x = Year + 1959, y = value, group = jitter, color = Hessian), lwd = 1) +
#'     geom_line(final_mod, mapping = aes(x = Year + 1959, y = value), color = "black", lwd = 1.3 , lty = 2) +
#'     facet_grid(Type~Region, scales = 'free',
#'                labeller = labeller(Region = function(x) paste0("Region ", x),
#'                                    Type = c("Recruitment" = "Age 2 Recruitment (millions)", "SSB" = 'SSB (kt)'))) +
#'     labs(x = "Year", y = "Value") +
#'     theme_bw(base_size = 20) +
#'     scale_color_manual(values = c("red", 'grey')) +
#'     geom_text(data = jit %>% filter(Type == 'SSB', Year == 1, jitter == 1),
#'               aes(x = Inf, y = Inf, label = paste("Proportion Converged: ", round(prop_converged$prop_conv, 3))),
#'               hjust = 1.1, vjust = 1.9, size = 6, color = "black")
#'
#'    # compare jitter of max gradient and hessian PD
#'    ggplot(jit, aes(x = jitter, y = jnLL, color = Max_Gradient, shape = Hessian)) +
#'      geom_point(size = 5, alpha = 0.3) +
#'      geom_hline(yintercept = min(sabie_rtmb_model$rep$jnLL), lty = 2, size = 2, color = "blue") +
#'      facet_wrap(~Hessian, labeller = labeller(
#'        Hessian = c("FALSE" = "non-PD Hessian", "TRUE" = 'PD Hessian')
#'      )) +
#'      scale_color_viridis_c() +
#'      theme_bw(base_size = 20) +
#'      theme(legend.position = "bottom") +
#'      guides(color = guide_colorbar(barwidth = 15, barheight = 0.5)) +
#'      labs(x = 'Jitter') +
#'      geom_text(data = jit %>% filter(Hessian == TRUE, Year == 1, jitter == 1),
#'                aes(x = Inf, y = Inf, label = paste("Proportion Converged: ", round(prop_converged$prop_conv, 3))),
#'                hjust = 1.1, vjust = 1.9, size = 6, color = "black")
#'}
do_jitter <- function(data,
                      parameters,
                      mapping,
                      random = NULL,
                      sd,
                      n_jitter,
                      n_newton_loops,
                      do_par,
                      n_cores
                      ) {

  jitter_all <- data.frame()

  obj <- RTMB::MakeADFun(cmb(SPoCK_rtmb, data),
                         parameters = parameters,
                         map = mapping,
                         random = random,
                         silent = TRUE)

  if(do_par == FALSE) {
    for(i in 1:n_jitter) {

      # jitter original parameters (additive normal draws)
      jitter_pars <- obj$par + stats::rnorm(length(obj$par), 0, sd)

      # Now, optimize the function
      optim <- stats::nlminb(jitter_pars,
                             obj$fn,
                             obj$gr,
                             control = list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15))

      # newton steps
      try_improve <- tryCatch(expr =
                                for(j in 1:n_newton_loops) {
                                  g = as.numeric(obj$gr(optim$par))
                                  h = optimHess(optim$par, fn = obj$fn, gr = obj$gr)
                                  optim$par = optim$par - solve(h,g)
                                  optim$objective = obj$fn(optim$par)
                                }
                              , error = function(e){e}, warning = function(w){w})

      obj$rep <- obj$report(obj$env$last.par.best) # Get report
      obj$sd_rep <- RTMB::sdreport(obj) # Get sd report

      # put jitter results into a dataframe
      jitter_ts_df <- reshape2::melt(obj$rep$SSB) %>%
        dplyr::rename(Region = Var1, Year = Var2) %>%
        dplyr::mutate(Type = 'SSB') %>%
        dplyr::bind_rows(reshape2::melt(obj$rep$Rec) %>%
                           dplyr::rename(Region = Var1, Year = Var2) %>%
                           dplyr::mutate(Type = 'Recruitment')) %>%
        dplyr::mutate(jitter = i,
                      Hessian = obj$sd_rep$pdHess,
                      jnLL = obj$rep$jnLL,
                      Max_Gradient = max(abs(obj$sd_rep$gradient.fixed)))

      jitter_all <- rbind(jitter_all, jitter_ts_df) # bind dataframes

    } # end i loop
  } # don't parrallelize

  if(do_par == TRUE) {

    plan(multisession, workers = n_cores) # set up cores

    with_progress({

      p <- progressor(along = 1:n_jitter) # progress bar

      jitter_all <- future_lapply(1:n_jitter, function(i) {

        # make obj
        obj <- RTMB::MakeADFun(cmb(SPoCK_rtmb, data), parameters = parameters,  map = mapping, random = random, silent = TRUE)

        # Jitter original parameters
        jitter_pars <- obj$par + stats::rnorm(length(obj$par), 0, sd)

        # Optimize function
        optim <- stats::nlminb(jitter_pars,
                               obj$fn,
                               obj$gr,
                               control = list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15))

        # Newton steps
        try_improve <- tryCatch({
          for (j in 1:n_newton_loops) {
            g <- as.numeric(obj$gr(optim$par))
            h <- optimHess(optim$par, fn = obj$fn, gr = obj$gr)
            optim$par <- optim$par - solve(h, g)
            optim$objective <- obj$fn(optim$par)
          }
        }, error = function(e) e, warning = function(w) w)

        # get reports
        obj$rep <- obj$report(obj$env$last.par.best)
        obj$sd_rep <- RTMB::sdreport(obj)

        # put jitter results into a dataframe
        jitter_ts_df <- reshape2::melt(obj$rep$SSB) %>%
          dplyr::rename(Region = Var1, Year = Var2) %>%
          dplyr::mutate(Type = 'SSB') %>%
          dplyr::bind_rows(reshape2::melt(obj$rep$Rec) %>%
                             dplyr::rename(Region = Var1, Year = Var2) %>%
                             dplyr::mutate(Type = 'Recruitment')) %>%
          dplyr::mutate(jitter = i,
                        Hessian = obj$sd_rep$pdHess,
                        jnLL = obj$rep$jnLL,
                        Max_Gradient = max(abs(obj$sd_rep$gradient.fixed)))

        p() # update progress

        jitter_ts_df

      }, future.seed = TRUE) %>% bind_rows() # bine rows to combine results

      plan(sequential)  # Reset

    })

  } # end dor par

  return(jitter_all)
}
