#' Get equilibrium age structure for simulation via n iterations
#'
#' @param n_regions Number of reginos
#' @param n_sexes Number of sexes
#' @param n_ages Number of ages
#' @param init_iter Number of initial iterations to do
#' @param NAA Initial numbers at age array dimensioned by year, region, age, sex, and 1
#' @param M Natural mortality array dimensioned by year, region, age, sex, and 1
#' @param init_F Initial fishing mortality array dimensioned by year, region, fleet, and 1
#' @param fish_sel Fishery selectivity array dimensioned by year, region, age, sex, fleet, and 1
#' @param r0 virgin or mean recruitment array dimensioned by year, region, and 1
#' @param rec_sexratio Recruitment sexratio array dimensioned by year, region, sex, and 1
#' @param movement_matrix Movement array, dimensioned by region from, region to, year, age, sex, and 1
#'
#' @keywords internal
#'
Equilibrium_PopSim_Iter <- function(n_regions = n_regions,
                                    n_sexes = n_sexes,
                                    n_ages = n_ages,
                                    init_iter = init_iter,
                                    NAA,
                                    M,
                                    init_F,
                                    fish_sel,
                                    r0,
                                    rec_sexratio,
                                    movement_matrix
                                    ) {

  NAA_nxt = NAA # set up variables here

  for(r in 1:n_regions) {
    for(s in 1:n_sexes) {
      tmp_cumsum_Z = cumsum(M[1,r,1:(n_ages-1),s,1] + init_F[1,r,1,1] * fish_sel[1,r,1:(n_ages-1),s,1,1]) # cumulative sum of Z with selectivity effects
      NAA[r,,s,1] = c(r0[1,r,1], r0[1,r,1] * exp(-tmp_cumsum_Z)) * rec_sexratio[1,r,s,1] # initial equilibrium age structure
    } # end s loop
  } # end r loop

  # Apply annual cycle and iterate to equilibrium
  for(i in 1:init_iter) {
    for(s in 1:n_sexes) {
      NAA_nxt[,1,s,1] = r0[1,,1] * rec_sexratio[1,,s,1] # recruitment

      # Movement
      if(do_recruits_move == 0) for(a in 2:n_ages) NAA[,a,s,1] = t(NAA[,a,s,1]) %*% movement_matrix[,,1,a,s,1] # Recruits don't move
      if(do_recruits_move == 1) for(a in 1:n_ages) NAA[,a,s,1] = t(NAA[,a,s,1]) %*% movement_matrix[,,1,a,s,1] # Recruits move

      # ageing and mortality
      NAA_nxt[,2:n_ages,s,1] = NAA[,1:(n_ages-1),s,1] * exp(-(M[1,,1:(n_ages-1),s,1] + (init_F[1,,,1] * fish_sel[1,,1:(n_ages-1),s,1,1])))
      # accumulate plus group
      NAA_nxt[,n_ages,s,1] = (NAA_nxt[,n_ages,s,1] * exp(-(M[1,,n_ages,s,1] + (init_F[1,,1,1] * fish_sel[1,,n_ages,s,1,1])))) +
                             (NAA[,n_ages,s,1] * exp(-(M[1,,n_ages,s,1] + (init_F[1,,1,1] * fish_sel[1,,n_ages,s,1,1]))))
      NAA = NAA_nxt # iterate to next cycle
    } # end s loop
  } # end i loop

  return(NAA)
} # end function
