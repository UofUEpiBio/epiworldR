#' Susceptible Exposed Infected Removed model (SEIR) with quarantine
#'
#' `ModelSEIRQuarantine` creates an SEIR model with detection, isolation,
#' contact tracing and quarantine. Its contacts come from its layers (see
#' [contact-layers]): a network, mixing by group, or both. Infected agents
#' who are detected isolate if willing, and their contacts in every layer are
#' traced and quarantined if willing.
#'
#' @param name String. Name of the virus
#' @param n Number of individuals in the population.
#' @param prevalence Double. Initial proportion of individuals with the virus.
#' @param transmission_prob Numeric scalar between 0 and 1. Probability that a
#' contact with an infected agent infects ("Transmission probability").
#' @param latent_period Numeric scalar. Mean days in Exposed: infected, not
#' yet infectious ("Latent period").
#' @param recovery_rate Numeric scalar between 0 and 1. Daily probability of
#' recovering ("Recovery rate").
#' @param hospitalization_rate Double. Daily probability of hospitalization
#' while infected ("Hospitalization rate").
#' @param hospitalization_period Double. Mean days in hospital
#' ("Hospitalization period").
#' @param detection_delay Double. Mean days infected before detection
#' ("Detection delay"; 0 detects at once, a negative value never).
#' @param quarantine_period Integer. Days in quarantine ("Quarantine period").
#' @param quarantine_willingness Double. Probability that an agent complies
#' with quarantine ("Quarantine willingness").
#' @param isolation_willingness Double. Probability that an agent isolates
#' when detected ("Isolation willingness").
#' @param isolation_period Integer. Days in isolation, from onset ("Isolation
#' period").
#' @param contact_tracing_prob Double. Probability that a contact is traced
#' ("Contact tracing probability").
#' @param contact_tracing_window Integer. Days before detection whose contacts
#' are traced ("Contact tracing window").
#' @export
#' @family Models
#' @concept general-models
#' @details
#' At the start of every run, the model turns contact tracing on for each of
#' its layers (see [enable_tracing()]). Only infected agents transmit;
#' isolated, quarantined and hospitalized agents are out of the mixing layers
#' altogether: they are not met, and their group's contacts go to the others.
#'
#' The quarantine and isolation processes can be turned off by specifying
#' `quarantine_period < 0` and `isolation_period < 0` respectively.
#'
#' The [initial_states] function allows the user to set the initial state of the
#' model. In particular, the user can specify how many of the non-infected
#' agents have been removed at the beginning of the simulation.
#' @returns
#' - The `ModelSEIRQuarantine` function returns a model of class [epiworld_model].
#' @aliases epiworld_seirquarantine
#'
#' @examples
#' N <- 9e3
#'
#' flu_model <- ModelSEIRQuarantine(
#'   name                   = "Flu",
#'   n                      = N,
#'   prevalence             = 1 / N,
#'   transmission_prob      = 0.1,
#'   latent_period          = 7,
#'   recovery_rate          = 1 / 7,
#'   hospitalization_rate   = 0.05,
#'   hospitalization_period = 7,
#'   detection_delay        = 3,
#'   quarantine_period      = 14,
#'   quarantine_willingness = 0.8,
#'   isolation_willingness  = 0.5,
#'   isolation_period       = 7,
#'   contact_tracing_prob   = 0.7,
#'   contact_tracing_window = 3
#' )
#'
#' # Three groups mixing by a contact matrix (20 contacts a day)
#' cmatrix <- matrix(
#'   c(0.9, 0.05, 0.05, 0.1, 0.8, 0.1, 0.1, 0.2, 0.7) * 20,
#'   byrow = TRUE,
#'   nrow  = 3
#' )
#' add_layer(
#'   flu_model,
#'   layer_mixing("community", groups = rep(0:2, each = 3e3), contact_matrix = cmatrix)
#' )
#'
#' run(flu_model, ndays = 100, seed = 331)
#' summary(flu_model)
#'
#' @seealso epiworld-methods
ModelSEIRQuarantine <- function(
  name,
  n,
  prevalence,
  transmission_prob,
  latent_period,
  recovery_rate,
  hospitalization_rate,
  hospitalization_period,
  detection_delay,
  quarantine_period,
  quarantine_willingness,
  isolation_willingness,
  isolation_period,
  contact_tracing_prob   = 1,
  contact_tracing_window = 4L
) {
  # Check input parameters
  stopifnot_string(name)
  stopifnot_int(n, lb = 1L)
  stopifnot_double(prevalence, lb = 0, ub = 1)
  stopifnot_double(transmission_prob, lb = 0, ub = 1)
  stopifnot_double(latent_period, lb = 0)
  stopifnot_double(recovery_rate, lb = 0, ub = 1)
  stopifnot_double(hospitalization_rate, lb = 0, ub = 1)
  stopifnot_double(hospitalization_period, lb = 0)
  stopifnot_double(detection_delay)
  stopifnot_int(quarantine_period)
  stopifnot_double(quarantine_willingness, lb = 0, ub = 1)
  stopifnot_double(isolation_willingness, lb = 0, ub = 1)
  stopifnot_int(isolation_period)
  stopifnot_double(contact_tracing_prob, lb = 0, ub = 1)
  stopifnot_int(contact_tracing_window, lb = 0)

  structure(
    ModelSEIRQuarantine_cpp(
      name, n, prevalence,
      transmission_prob, latent_period, recovery_rate,
      hospitalization_rate, hospitalization_period,
      detection_delay, quarantine_period,
      quarantine_willingness, isolation_willingness, isolation_period,
      contact_tracing_prob, contact_tracing_window
    ),
    class = c("epiworld_seirquarantine", "epiworld_model")
  )

}
