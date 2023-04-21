#' Virus
#' @param name of the virus
#' @param post_immunity,prob_infecting,prob_recovery,prob_death to be documented
#' @param virus_pos Positive integer. Index of the virus's position in the model.
#' @examples 
#' mseirconn <- ModelSEIRCONN(
#'   name                = "COVID-19",
#'   prevalence          = 0.01, 
#'   n                   = 10000,
#'   contact_rate        = 4, 
#'   incubation_days     = 7, 
#'   prob_transmission   = 0.5,
#'   prob_recovery       = 0.99
#' )
#' 
#' delta <- virus("Delta Variant", 0, .5, .2, .01)
#' virus_set_state(delta, 1, 3, 3)
#' 
#' add_virus(mseirconn, delta, .3)
#' 
#' run(mseirconn, ndays = 100, seed = 992)
#' mseirconn
#' @export
#' 
virus <- function(
    name,
    prob_infecting,
    prob_recovery = 0.5,
    prob_death    = 0.0,
    post_immunity = 1.0
    ) {
  
  structure(
    virus_cpp(
      name,
      prob_infecting,
      prob_recovery,
      prob_death,
      post_immunity
      ),
    class = "epiworld_virus"
  )
  
}

#' @export
#' @rdname virus
#' @param m An object of class `epiworld-model`.
#' @param v An object of class `epiworld_virus`
#' @param prevalence In the case of `add_virus`, a proportion, otherwise, an integer.
add_virus <- function(m, v, prevalence) UseMethod("add_virus")

#' @export
add_virus.epiworld_model <- function(m, v, prevalence) {
  
  add_virus_cpp(m, v, prevalence)
  invisible(m)
  
}

#' @export
add_virus.epiworld_sir <- function(m, v, prevalence) {
  
  virus_set_state(v, init = 1, end = 2, removed = 2)
  add_virus_cpp(m, v, prevalence)
  invisible(m)
  
}

#' @export
add_virus.epiworld_sirconn <- function(m, v, prevalence) {
  
  add_virus.epiworld_sir(m, v, prevalence)
  
}

#' @export
add_virus.epiworld_seir <- function(m, v, prevalence) {
  
  add_virus.epiworld_sir(m, v, prevalence)
  
}

#' @export
add_virus.epiworld_seirconn <- function(m, v, prevalence) {
  
  add_virus.epiworld_sir(m, v, prevalence)
  
}

#' @export
#' @rdname virus
add_virus_n <- function(m, v, prevalence) UseMethod("add_virus_n")

#' @export
add_virus_n.epiworld_model <- function(m, v, prevalence) {
  
  add_virus_n_cpp(m, v, prevalence)
  
  invisible(m)
}

#' @export
#' @rdname virus
#' @param init,end,removed states after acquiring a virus, removing a virus,
#' and removing the agent as a result of the virus, respectively.
virus_set_state <- function(v, init, end, removed) {
  
  if (!inherits(v, "epiworld_virus"))
    stop("-v- must be of class epiworld_virus.")
  
  virus_set_state_cpp(v, init, end, removed)
  invisible(v)
  
}

#' @export
#' @rdname virus
rm_virus <- function(m, virus_pos) {
  invisible(rm_virus_cpp(m, virus_pos))
}


#' @export
print.epiworld_virus <- function(x, ...) {
  print_virus_cpp(x)
}

