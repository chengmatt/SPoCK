#' GGPLOT theme for sablefish
#'
#' @return ggplot theme
#' @export theme_sablefish
#' @import ggplot2
theme_sablefish <- function() {
   theme_bw() +
    theme(legend.position = "top",
          strip.text = element_text(size = 17),
          title = element_text(size = 21, color = 'black'),
          axis.text = element_text(size = 15, color = "black"),
          axis.= element_text(size = 17, color = 'black'),
          legend.text = element_text(size = 15, color = "black"),
          legend.= element_text(size = 17, color = 'black'))
}


#' Function to fill in an n x n correlation AR(1) matrix
#'
#' @param n Number of bins
#' @param rho correaltion parameter
#'
#' @return correlation matrix for an ar1 process
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' mat <- get_AR1_CorrMat(10, 0.5)
#' }
get_AR1_CorrMat <- function(n, rho) {
  corrMatrix <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      # Calculate the correlation based on the lag distance
      corrMatrix[i, j] <- rho^(abs(i - j))
    } # end i
  } # end j
  return(corrMatrix)
}

#' Constant correlation matrix
#'
#' @param n Number of bins
#' @param rho correaltion parameter
#'
#' @return constant correlation matrix
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' mat <- get_Constant_CorrMat(10, 0.5)
#' }
get_Constant_CorrMat <- function(n, rho) {
  corrMatrix <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      if(i != j) corrMatrix[i, j] <- rho
      else corrMatrix[i, j] <- 1
    } # end i
  } # end j
  return(corrMatrix)
}

#' For combining a parameter and data list in RTMB so a data object can be explicitly defined
#'
#' @param f Parameter list
#' @param d Data list
#' @keywords internal
#' @examples
#' \dontrun{
#'   obj <- RTMB::MakeADFun(cmb(sabie_RTMB, data), parameters = parameters, map = mapping, random = random, silent = TRUE)
#' }
cmb <- function(f, d) {
  function(p) f(p, d)
}
