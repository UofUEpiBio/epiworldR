#' Virus design
#' 
#' Viruses can be considered to be anything that can be transmitted (e.g.,
#' diseases, as well as ideas.) Most models in epiworldR can feature multiple 
#' viruses.
#' 
#' @param name of the virus
#' @param post_immunity Numeric scalar. Post immunity (prob of re-infection).
#' @param prob_infecting Numeric scalar. Probability of infection (transmission).
#' @param prob_recovery Numeric scalar. Probability of recovery.
#' @param prob_death Numeric scalar. Probability of death.
#' @param virus_pos Positive integer. Index of the virus's position in the model.
#' @details
#' The [virus()] function can be used to initialize a virus. Virus features can
#' then be modified using the functions `set_prob_*`.
#' 
#' The function [virus_fun_logit()] creates a "virus function" that can be
#' evaluated for transmission, recovery, and death. As the name sugests, it
#' computes those probabilities using a logit function (see examples).
#' 
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
    post_immunity = -1.0
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
print.epiworld_virus <- function(x, ...) {
  print_virus_cpp(x)
}

stopifnot_virus <- function(virus) {
  if (!inherits(virus, "epiworld_virus")) {
    stop(
      "The -virus- object must be of class \"epiworld_virus\". ",
      "The object passed to the function is of class(es): ", 
      paste(class(virus), collapse = ", ")
    )
  }
}

stopifnot_vfun <- function(vfun) {
  if (!inherits(vfun, "epiworld_virus_fun")) {
    stop(
      "The -vfun- object must be of class \"epiworld_virus_fun\". ",
      "The object passed to the function is of class(es): ", 
      paste(class(vfun), collapse = ", ")
    )
  }
}

#' @export
#' @details
#' The name of the `epiworld_virus` object can be manipulated with the functions
#' [set_name_virus()] and [get_name_virus()].
#' 
#' @rdname virus
set_name_virus <- function(virus, name) {
  stopifnot_virus(virus)
  invisible(set_name_virus_cpp(virus, name))
}

#' @export
#' @rdname virus
get_name_virus <- function(virus) {
  stopifnot_virus(virus)
  get_name_virus(virus)
}
  
# Virus add --------------------------------------------------------------------

#' @export
#' @rdname virus
#' @param model An object of class `epiworld-model`.
#' @param virus An object of class `epiworld_virus`
#' @param prevalence In the case of `add_virus`, a proportion, otherwise, an integer.
add_virus <- function(model, virus, prevalence) UseMethod("add_virus")

#' @export
add_virus.epiworld_model <- function(model, virus, prevalence) {
  
  stopifnot_virus(virus)
  
  add_virus_cpp(model, virus, prevalence)
  invisible(model)
  
}

#' @export
add_virus.epiworld_sir <- function(model, virus, prevalence) {
  
  stopifnot_virus(virus)
  virus_set_state(virus, init = 1, end = 2, removed = 2)
  invisible(add_virus_cpp(model, virus, prevalence))
  
}

#' @export
add_virus.epiworld_sirconn <- function(model, virus, prevalence) {
  
  stopifnot_virus(virus)
  add_virus.epiworld_sir(model, virus, prevalence)
  
}

#' @export
add_virus.epiworld_seir <- function(model, virus, prevalence) {
  
  stopifnot_virus(virus)
  virus_set_state(virus, init = 1, end = 3, removed = 3)
  invisible(add_virus_cpp(model, virus, prevalence))
  
}

#' @export
add_virus.epiworld_seirconn <- function(model, virus, prevalence) {
  
  stopifnot_virus(virus)
  add_virus.epiworld_seir(model, virus, prevalence)
  
}

#' @export
#' @rdname virus
add_virus_n <- function(model, virus, prevalence) UseMethod("add_virus_n")

#' @export
add_virus_n.epiworld_model <- function(model, virus, prevalence) {
  
  stopifnot_virus(virus)
  invisible(add_virus_n_cpp(model, virus, prevalence))
  
}

#' @export
add_virus_n.epiworld_sir <- function(model, virus, prevalence) {
  
  stopifnot_virus(virus)
  virus_set_state(model, init = 1, end = 2, removed = 2)
  invisible(add_virus_n_cpp(model, virus, prevalence))
  
}

#' @export
add_virus_n.epiworld_sirconn <- function(model, virus, prevalence) {
  
  stopifnot_virus(virus)
  add_virus_n.epiworld_sir(model, virus, prevalence)
  
}

#' @export
add_virus_n.epiworld_seir <- function(model, virus, prevalence) {
  
  stopifnot_virus(virus)
  virus_set_state(model, init = 1, end = 3, removed = 3)
  invisible(add_virus_n_cpp(model, virus, prevalence))
  
}

#' @export
add_virus_n.epiworld_seirconn <- function(model, virus, prevalence) {
  
  stopifnot_virus(virus)
  add_virus_n.epiworld_seir(model, virus, prevalence)
  
}

# Virus MISC -------------------------------------------------------------------

#' @export
#' @rdname virus
#' @param init,end,removed states after acquiring a virus, removing a virus,
#' and removing the agent as a result of the virus, respectively.
virus_set_state <- function(virus, init, end, removed) {
  
  stopifnot_virus(virus)
  invisible(virus_set_state_cpp(virus, init, end, removed))
  
}

#' @export
#' @rdname virus
rm_virus <- function(model, virus_pos) {
  
  stopifnot_model(model)
  invisible(rm_virus_cpp(model, virus_pos))
  
}

# Virus functions --------------------------------------------------------------

#' @export
#' @param vars Integer vector. Indices (starting from 0) of the positions of the
#' variables used to compute the logit probability.
#' @param coefs Numeric vector. Of the same length of `vars`, is a vector of
#' coefficients associated to the logit probability.
#' @rdname virus
#' @examples
#' # Using the logit function --------------
#' sir <- ModelSIR(
#'   name = "COVID-19", prevalence = 0.01, 
#'   infectiousness = 0.9, recovery = 0.1
#'   )
#' 
#' # Adding a small world population
#' agents_smallworld(
#'   sir,
#'   n = 10000,
#'   k = 5,
#'   d = FALSE,
#'   p = .01
#' )
#' 
#' run(sir, ndays = 50, seed = 11)
#' plot(sir)
#' 
#' # And adding features
#' dat <- cbind(
#'   female = sample.int(2, 10000, replace = TRUE) - 1,
#'   x      = rnorm(10000)
#' )
#' 
#' set_agents_data(sir, dat)
#' 
#' # Creating the logit function
#' vfun <- virus_fun_logit(
#'   vars  = c(0L, 1L),
#'   coefs = c(-1, 1),
#'   model = sir
#' )
#' 
#' # The infection prob is lower
#' hist(plogis(dat %*% rbind(-1,1)))
#' 
#' vfun # printing
#' 
#' set_prob_infecting_fun(
#'   virus = get_virus(sir, 0),
#'   model = sir,
#'   vfun  = vfun
#'   )
#'   
#' run(sir, ndays = 50, seed = 11)
#' plot(sir)
#' 
#' 
virus_fun_logit <- function(vars, coefs, model) {
  
  stopifnot_model(model)
  
  structure(
    virus_fun_logit_cpp(vars, coefs, model),
    class = "epiworld_virus_fun",
    builder = "virus_fun_logit",
    vars    = vars,
    coefs   = coefs,
    model   = model
  )
  
}

#' @export
print.epiworld_virus_fun <- function(x, ...) {
  
  cat("An epiworld_virus_function object.\n")
  cat("(model: ", get_name(attr(x, "model")), ")\n", sep = "")
  cat("This function was built using -virus_fun_logit()-. and it features ")
  cat("the following coefficients:\n")
  cat(
    paste(sprintf(
      " % 2i: %5.2f",
      attr(x, "vars"),
      attr(x, "coefs")
      ), collapse = "\n"
    ), "\n"
  )
  
  invisible(x)
  
}


#' @export
#' @param prob Numeric scalar. A probability (between zero and one).
#' @rdname virus
set_prob_infecting <- function(virus, prob) {
  
  stopifnot_virus(virus)
  set_prob_infecting_cpp(virus, prob)
  
}

#' @export
#' @param param Character scalar. Name of the parameter featured in `model` that
#' will be added to the virus (see details.)
#' @details
#' In the case of `set_prob_infecting_ptr`, `set_prob_recovery_ptr`, and
#' `set_prob_death_ptr`, the corresponding parameters is passed as a pointer to
#' the virus. The implication of using pointers is that the values will be
#' read directly from the `model` object, so changes will be reflected.
#' 
#' @rdname virus
set_prob_infecting_ptr <- function(virus, model, param) {
  
  stopifnot_virus(virus)
  stopifnot_model(model)
  invisible(set_prob_infecting_ptr_cpp(virus, model, param))
  
}

#' @export
#' @param vfun An object of class `epiworld_virus_fun`.
#' @rdname virus
set_prob_infecting_fun <- function(virus, model, vfun) {
  
  stopifnot_virus(virus)
  stopifnot_model(model)
  stopifnot_vfun(vfun)
  invisible(set_prob_infecting_fun_cpp(virus, model, vfun))
  
}

#' @export
#' @rdname virus
set_prob_recovery <- function(virus, prob) {
  
  stopifnot_virus(virus)
  set_prob_recovery_cpp(virus, prob)
  
}

#' @export
#' @rdname virus
set_prob_recovery_ptr <- function(virus, model, param) {
  
  stopifnot_virus(virus)
  stopifnot_model(model)
  set_prob_recovery_ptr_cpp(virus, model, param)
  
}

#' @export
#' @rdname virus
set_prob_recovery_fun <- function(virus, model, vfun) {
  
  stopifnot_virus(virus)
  stopifnot_model(model)
  stopifnot_vfun(vfun)
  set_prob_recovery_fun_cpp(virus, model, vfun)
}

#' @export
#' @rdname virus
set_prob_death <- function(virus, prob) {
  
  stopifnot_virus(virus)
  set_prob_death_cpp(virus, prob)
  
}

#' @export
#' @rdname virus
set_prob_death_ptr <- function(virus, model, param) {
  
  stopifnot_virus(virus)
  stopifnot_model(model)
  set_prob_death_ptr_cpp(virus, model, param)
  
}

#' @export
#' @rdname virus
set_prob_death_fun <- function(virus, model, vfun) {
  
  stopifnot_virus(virus)
  stopifnot_model(model)
  stopifnot_vfun(vfun)
  set_prob_death_fun_cpp(virus, model, vfun)
  
}


