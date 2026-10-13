# Shared recipe for the six standard network models. Public constructors keep
# their argument validation and S3 classes; only assembly happens in R.
network_model_builder <- function(
  kind, name, prevalence, transmission_rate, recovery_rate,
  incubation_days = NULL, death_rate = NULL
) {
  model <- Model()
  latent <- kind %in% c("SEIR", "SEIRD")
  reinfection <- kind %in% c("SIS", "SISD")
  mortality <- !is.null(death_rate)

  add_state(model, "Susceptible", update_fun_susceptible())
  if (latent)
    add_state(model, "Exposed", update_fun_incubation(2L))
  add_state(model, "Infected", if (kind == "SEIR") {
    update_fun_recovery("Recovery rate")
  } else {
    update_fun_infected()
  })
  if (!reinfection)
    add_state(model, if (latent) "Removed" else "Recovered")
  if (mortality)
    add_state(model, "Deceased")

  # Preserve parameter order as well as names and values.
  if (kind %in% c("SIR", "SIRD")) {
    add_param(model, "Recovery rate", recovery_rate)
    add_param(model, "Transmission rate", transmission_rate)
  } else {
    add_param(model, "Transmission rate", transmission_rate)
    if (latent)
      add_param(model, "Incubation days", incubation_days)
    add_param(model, "Recovery rate", recovery_rate)
  }
  if (mortality)
    add_param(model, "Death rate", death_rate)

  pathogen <- virus(name, prevalence, TRUE, transmission_rate, recovery_rate)
  acquired <- 1L
  recovered <- if (reinfection) 0L else if (latent) 3L else 2L
  deceased <- if (!mortality) {
    recovered
  } else if (latent) {
    4L
  } else if (reinfection) {
    2L
  } else {
    3L
  }
  virus_set_state(pathogen, acquired, recovered, deceased)
  set_prob_infecting_ptr(pathogen, model, "Transmission rate")
  set_prob_recovery_ptr(pathogen, model, "Recovery rate")
  if (latent)
    set_incubation_ptr(pathogen, model, "Incubation days")
  if (kind == "SISD") {
    # Preserve the existing C++ model's fixed death probability. Correcting
    # its ignored Death rate parameter is a separate behavior change.
    set_prob_death(pathogen, 0.01)
  } else if (mortality) {
    set_prob_death_ptr(pathogen, model, "Death rate")
  }
  add_virus(model, pathogen)

  labels <- c("Susceptible", if (latent) "Exposed", "Infected",
    if (reinfection) "Susceptible" else if (latent) "Removed" else "Recovered",
    if (mortality) "Deceased")
  set_name(model, paste0(paste(labels, collapse = "-"), " (", kind, ")"))
  if (!reinfection)
    set_initial_states_recipe_cpp(model, kind)
  model
}
