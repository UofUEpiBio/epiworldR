#' SIR Logistic model 
#' @param vname 
#' @param data 
#' @param coefs_infect 
#' @param coefs_recover 
#' @param coef_infect_cols 
#' @param coef_recover_cols 
#' @param prevalence 
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

