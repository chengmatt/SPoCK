#' Constructor algorithin for correlations within ages, years, and cohort (Work in Progress)
#'
#' @param n_ages Number of ages
#' @param n_yrs Number of years
#' @param pcorr_age correlations for age
#' @param pcorr_year correaltions for year
#' @param pcorr_cohort correlaitons for cohort
#' @param ln_var_value log space variance
#' @param Var_Type variance type == 0, marginal (stationary and slower run time), == 1 conditional (non-statationary, faster run time)
#'
#' @importFrom Matrix sparseMatrix
#' @returns Sparse precision matrix dimensioned by n_ages * n_years, n_ages * n_years
#' @export Get_3d_precision
#'
Get_3d_precision <- function(n_ages,n_yrs, pcorr_age, pcorr_year, pcorr_cohort, ln_var_value, Var_Type){

    "c" <- RTMB::ADoverload("c")
    "[<-" <- RTMB::ADoverload("[<-")

    index = expand.grid(seq_len(n_ages), seq_len(n_yrs)) # create index combinations to loop through
    i = j = x = numeric(0) # initialize posiiton to fill in precision matrix
    var_value = exp(ln_var_value) # transform to normal space

    for(n in 1:nrow(index)){
      age = index[n,1] # get age index out of all index combinations
      year = index[n,2] # get year index out of all index combinations
      if(age > 1 ){
        i = c(i, n)
        j = c(j, which(index[,1] == (age-1) & index[,2] == year))
        x = c(x, pcorr_year) # year correaltion indexing
      }
      if(year > 1){
        i = c(i, n)
        j = c(j, which(index[,1]==age & index[,2]==(year-1)) )
        x = c(x, pcorr_age) # age correlation indexing
      }
      if( age>1 & year>1 ){
        i = c(i, n)
        j = c(j, which(index[,1]==(age-1) & index[,2] == (year-1)) )
        x = c(x, pcorr_cohort) # cohort correlation indexing
      }
    } # end n loop

    # create B path matrix
    B = matrix(0, nrow = n_ages * n_yrs, ncol = n_ages * n_yrs)
    B[cbind(i, j)] = x
    B = as(B, "sparseMatrix")

    # identity matrix
    I = as(diag(1, n_ages * n_yrs, n_ages * n_yrs), "sparseMatrix")

    # Solve Omega recursively for stationary variance (accumulator function)
    if(Var_Type == 0) {
      L = solve(I-B) # solve to get accumulator function for stationary variance
      d = rep(0, nrow(index))
      for(n in 1:nrow(index) ){
        if(n==1){
          d[n] = var_value
        }else{
          cumvar = sum(L[n,seq_len(n-1)] * d[seq_len(n-1)] * L[n,seq_len(n-1)])
          d[n] = (var_value-cumvar) / L[n,n]^2
        }
      } # end n loop
    } # end marginal variance (stationary variance)

    if(Var_Type == 1) d = var_value # conditional variance (non-stationary variance)

    # omega matrix
    Omega_inv = diag(1/d, n_ages * n_yrs, n_ages * n_yrs)
    Q = as((I-t(B)) %*% Omega_inv %*% (I-B), "sparseMatrix") # solve for precision

    return(Q)
  }



