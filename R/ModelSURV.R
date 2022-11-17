#' Susceptible Infected Susceptible model (SURV)
#'
#' @param name String. Name of the virus
#' @param prevalence Initial number of individuals with the virus.
#' @param efficacy_vax Double. Efficacy of the vaccine 
#' (1 - P(acquire the disease)).
#' @param latent_period Double. Shape parameter of a 'Gamma(latent_period, 1)'
#' distribution. This coincides with the expected number of latent days.
#' @param infect_period Double. Shape parameter of a 'Gamma(infected_period, 1)'
#' distribution. This coincides with the expected number of infectious days.
#' @param prob_symptoms Double. Probability of generating symptoms.
#' @param prop_vaccinated Double. Probability of vaccination. Coincides with
#' the initial prevalence of vaccinated individuals.
#' @param prop_vax_redux_transm Double. Factor by which the vaccine reduces
#' transmissibility.
#' @param prop_vax_redux_infect Double. Factor by which the vaccine reduces
#' the chances of becoming infected.
#' @param surveillance_prob Double. Probability of testing an agent.
#' @param prob_transmission Double. Raw transmission probability.
#' @param prob_death Double. Raw probability of death for symptomatic 
#' individuals.
#' @param prob_noreinfect Double. Probability of no re-infection.
#' @param x Object of class SIS. 
#' @param ... Currently ignore. 
#' @param n Number of individuals in the population.
#' @param k Number of ties in the small world network.
#' @param d Logical scalar. Whether the graph is directed or not.
#' @param p Probability of rewiring.
#' @export
#' @family Models
#' @aliases epiworld_surv
#' @examples
#' model_surv <- ModelSURV(
#'   name                  = "COVID-19",
#'   prevalence            = 20,
#'   efficacy_vax          = 0.6,
#'   latent_period         = 4,
#'   infect_period         = 5,
#'   prob_symptoms         = 0.5,
#'   prop_vaccinated       = 0.7,
#'   prop_vax_redux_transm = 0.8,
#'   prop_vax_redux_infect = 0.95,
#'   surveillance_prob     = 0.1,
#'   prob_transmission     = 0.2,
#'   prob_death            = 0.001,
#'   prob_noreinfect       = 0.5
#' )
#'
#' # Adding a small world population
#' agents_smallworld(
#'   model_surv,
#'   n = 10000,
#'   k = 5,
#'   d = FALSE,
#'   p = .01
#'   )
#'   
#' # Initializing
#' init(model_surv, days = 100, seed = 1912)
#' # Running and printing
#' run(model_surv)
#' model_surv 
#' 

ModelSURV <- function(
    name, prevalence, efficacy_vax, latent_period, infect_period, prob_symptoms, 
    prop_vaccinated, prop_vax_redux_transm, prop_vax_redux_infect, 
    surveillance_prob, prob_transmission, prob_death, prob_noreinfect  
) {
  
  structure(
    ModelSURV_cpp(name, prevalence, efficacy_vax, latent_period, infect_period, 
                  prob_symptoms, prop_vaccinated, prop_vax_redux_transm, 
                  prop_vax_redux_infect, surveillance_prob, prob_transmission, 
                  prob_death, prob_noreinfect),
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
  print_surv(x)
}

#' @rdname ModelSURV
#' @export
agents_smallworld.epiworld_surv <- function(m, n, k, d, p) {
  agents_smallworld_surv(m, n, k, d, p)
}

#' @rdname ModelSURV
#' @export
run.epiworld_surv <- function(m) {
  run_surv(m)
}

#' @rdname ModelSURV
#' @export
plot.epiworld_surv <- function(x, ...) { # col = NULL
 plot_epi(x, main = "SURV Model", ...)
}
