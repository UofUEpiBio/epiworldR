#' SIR Logistic model
#' @param vname Name of the virus.
#' @param data A numeric matrix with `n` rows.
#' @param coefs_infect Numeric vector. Coefficients associated to infect.
#' @param coefs_recover Numeric vector. Coefficients associated to recover.
#' @param coef_infect_cols Integer vector. Columns in the coeficient.
#' @param coef_recover_cols Integer vector. Columns in the coeficient.
#' @param prob_infection Numeric scalar. Baseline probability of infection.
#' @param recovery_rate  Numeric scalar. Baseline probability of recovery.
#' @param prevalence Numeric scalar. Prevalence (initial state) in proportion.
#' @section Model diagram:
#' ![](diagrams/sirlogit.png "SIR Logit Diagram")
#' @export
#' @concept general-models
#' @returns
#' - The `ModelSIRLogit` function returns a model of class [epiworld_model].
#' @examples
#' set.seed(2223)
#' n <- 100000
#'
#' # Creating the data to use for the "ModelSIRLogit" function. It contains
#' # information on the sex of each agent and will be used to determine
#' # differences in disease progression between males and females. Note that
#' # the number of rows in these data are identical to n (100000).
#' X <- cbind(
#'   Intercept = 1,
#'   Female    = sample.int(2, n, replace = TRUE) - 1
#' )
#'
#' # Declare coefficients for each sex regarding transmission_rate and recovery.
#' coef_infect  <- c(.1, -2, 2)
#' coef_recover <- rnorm(2)
#'
#' # Feed all above information into the "ModelSIRLogit" function.
#' model_logit <- ModelSIRLogit(
#'   "covid2",
#'   data = X,
#'   coefs_infect      = coef_infect,
#'   coefs_recover     = coef_recover,
#'   coef_infect_cols  = 1L:ncol(X),
#'   coef_recover_cols = 1L:ncol(X),
#'   prob_infection = .8,
#'   recovery_rate = .3,
#'   prevalence = .01
#' )
#'
#' agents_smallworld(model_logit, n, 8, FALSE, .01)
#'
#' run(model_logit, 50)
#'
#' plot(model_logit)
#'
#' # Females are supposed to be more likely to become infected.
#' rn <- get_reproductive_number(model_logit)
#'
#' # Probability of infection for males and females.
#' (table(
#'   X[, "Female"],
#'   (1:n %in% rn$source)
#' ) |> prop.table())[, 2]
#'
#' # Looking into the individual agents.
#' get_agents(model_logit)

#' @family Models
ModelSIRLogit <- function(
  vname,
  data,
  coefs_infect,
  coefs_recover,
  coef_infect_cols,
  coef_recover_cols,
  prob_infection,
  recovery_rate,
  prevalence
) {
  # Check input parameters
  stopifnot_string(vname)
  stopifany_na(data)
  stopifnot_numvector(coefs_infect)
  stopifnot_numvector(coefs_recover)
  stopifnot_numvector(coef_infect_cols)
  stopifnot_numvector(coef_recover_cols)
  stopifnot_double(prob_infection)
  stopifnot_double(recovery_rate)
  stopifnot_double(prevalence)

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
      recovery_rate,
      prevalence
    ),
    class = c("epiworld_sir", "epiworld_model")
  )

}
