#' SIR Logistic model 
#' @param vname Name of the virus
#' @param data A numeric matrix with `n` rows.
#' @param coefs_infect Numeric vector. Coefficients associated to infect.
#' @param coefs_recover Numeric vector. Coefficients associated to recover.
#' @param coef_infect_cols Integer vector. Columns in the coeficient.
#' @param coef_recover_cols Integer vector. Columns in the coeficient.
#' @param prob_infection Numeric scalar. Baseline probability of infection.
#' @param prob_recovery  Numeric scalar. Baseline probability of recovery.
#' @param prevalence Numeric scalar. Prevalence (initial state) in proportion.
#'
#' @export
#' @family Models
ModelSIRLogit <- function(
  vname,
  data,
  coefs_infect,
  coefs_recover,
  coef_infect_cols,
  coef_recover_cols,
  prob_infection,
  prob_recovery,
  prevalence
) {
  
  structure(
    ModelSIRLogit_cpp(
      vname,
      data,
      ncol(data),
      coefs_infect,
      coefs_recover,
      coef_infect_cols - 1L,
      coef_recover_cols - 1L,
      prob_infection,
      prob_recovery,
      prevalence
    ),
    class = c("epiworld_sir", "epiworld_model")
  )
  
}

