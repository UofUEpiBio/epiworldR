#' Susceptible Infected Susceptible model (SURV)
#'
#' @param name Name of the virus
#' @param efficacy_vax to be explained
#' @param latent_period to be explained
#' @param infect_period to be explained
#' @param prob_symptoms to be explained
#' @param prop_vaccinated to be explained
#' @param prop_vax_redux_transm to be explained
#' @param prop_vax_redux_infect to be explained
#' @param surveillance_prob to be explained
#' @param prob_transmission to be explained
#' @param prob_death to be explained
#' @param prob_noreinfect to be explained
#' @param prevalence a number
#' @param m,days,seed,x,...,n,k,d,p to be explained
#' @export
#' @family Models
#' @aliases epiworld_surv
ModelSURV <- function(
    name, prevalence, efficacy_vax, latent_period, infect_period, prob_symptoms, prop_vaccinated, prop_vax_redux_transm, 
    prop_vax_redux_infect, surveillance_prob, prob_transmission, prob_death, prob_noreinfect  
) {
  
  structure(
    ModelSURV_cpp(name, prevalence, efficacy_vax, latent_period, infect_period, prob_symptoms, prop_vaccinated, prop_vax_redux_transm, 
    prop_vax_redux_infect, surveillance_prob, prob_transmission, prob_death, prob_noreinfect),
    class = "epiworld_surv"
  )
  
}

#' @rdname ModelSURV
#' @export
init.epiworld_surv <- function(m, days, seed) {
  init_surv(m, days, seed)
}

#' @rdname ModelSURV
#' @export
print.epiworld_surv <- function(x, ...) {
  print_sir(x)
}

#' @rdname ModelSURV
#' @export
agents_smallworld.epiworld_surv <- function(m, n, k, d, p) {
  agents_smallworld_sir(m, n, k, d, p)
}

#' @rdname ModelSURV
#' @export
run.epiworld_surv <- function(m) {
  run_sir(m)
}

#' @rdname ModelSURV
#' @export
plot.epiworld_surv <- function(x, ...) { # col = NULL
 plot_epi(x, main = "SURV Model", ...)
}
