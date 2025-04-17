#' Run RTMB model
#'
#' @param data Data list
#' @param parameters Parameter list
#' @param mapping Mapping list
#' @param random Character of random effects to integrate out
#' @param newton_loops Number of newton loops to run to get gradients down
#'
#' @return Returns a list object that is optimized, with results outputted from the RTMB model
#' @export fit_model
#'
#' @examples
#' \dontrun{
#'model <- fit_model(data,
#'                  parameters,
#'                  mapping,
#'                  random = NULL,
#'                  newton_loops = 3)
#'}
fit_model <- function(data,
                      parameters,
                      mapping,
                      random = NULL,
                      newton_loops = 3
                      ) {

  # make AD model function
  obj <- RTMB::MakeADFun(cmb(SPoCK_rtmb, data), parameters = parameters, map = mapping, random = random)

  # Now, optimize the function
  optim <- stats::nlminb(obj$par, obj$fn, obj$gr,
                               control = list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15))
  # newton steps
  try_improve <- tryCatch(expr =
                            for(i in 1:newton_loops) {
                              g = as.numeric(obj$gr(optim$par))
                              h = optimHess(optim$par, fn = obj$fn, gr = obj$gr)
                              optim$par = optim$par - solve(h,g)
                              optim$objective = obj$fn(optim$par)
                            }
                          , error = function(e){e}, warning = function(w){w})

  # save report
  obj$rep <- obj$report(obj$env$last.par.best)

  return(obj)
}
