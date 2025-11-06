#' Virus design
#'
#' Viruses can be considered to be anything that can be transmitted (e.g.,
#' diseases, as well as ideas.) Most models in epiworldR can feature multiple
#' viruses.
#'
#' @param name of the virus
#' @param post_immunity Numeric scalar. Post immunity (prob of re-infection).
#' @param prob_infecting Numeric scalar. Probability of infection (transmission).
#' @param recovery_rate Numeric scalar. Probability of recovery.
#' @param prob_death Numeric scalar. Probability of death.
#' @param virus_pos Positive integer. Index of the virus's position in the model.
#' @param incubation Numeric scalar. Incubation period (in days) of the virus.
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
#'   transmission_rate   = 0.5,
#'   recovery_rate       = 0.99
#' )
#'
#' delta <- virus(
#'   "Delta Variant", 0, .5, .2, .01, prevalence = 0.3, as_proportion = TRUE
#' )
#'
#' # Adding virus and setting/getting virus name
#' add_virus(mseirconn, delta)
#' set_name_virus(delta, "COVID-19 Strain")
#' get_name_virus(delta)
#'
#' run(mseirconn, ndays = 100, seed = 992)
#' mseirconn
#'
#' rm_virus(mseirconn, 0) # Removing the first virus from the model object
#' set_distribution_virus(delta, distribute_virus_randomly(100, as_proportion = FALSE))
#' add_virus(mseirconn, delta)
#'
#' # Setting parameters for the delta virus manually
#' set_prob_infecting(delta, 0.5)
#' set_prob_recovery(delta, 0.9)
#' set_prob_death(delta, 0.01)
#' run(mseirconn, ndays = 100, seed = 992) # Run the model to observe changes
#'
#' # If the states were (for example):
#' # 1: Infected
#' # 2: Recovered
#' # 3: Dead
#' delta2 <- virus(
#'   "Delta Variant 2", 0, .5, .2, .01, prevalence = 0, as_proportion = TRUE
#' )
#' virus_set_state(delta2, 1, 2, 3)
#' @export
#' @aliases epiworld_virus
virus <- function(
  name,
  prevalence,
  as_proportion,
  prob_infecting,
  recovery_rate = 0.5,
  prob_death    = 0.0,
  post_immunity = -1.0,
  incubation    = 7.0
) {

  uses_deprecated <- FALSE
  if (missing(prevalence)) {

    warning(
      "Starting version 0.3-0, the 'prevalence' argument is required.",
      " It will be set to be 0.5. Next versions will fail with an error."
    )

    prevalence <- 0.5
    as_proportion <- TRUE
    uses_deprecated <- TRUE

  }

  structure(
    virus_cpp(
      name,
      prevalence,
      as_proportion,
      prob_infecting,
      recovery_rate,
      prob_death,
      post_immunity,
      incubation
    ),
    class = "epiworld_virus",
    uses_deprecated = uses_deprecated,
    deprecated_args = list(
      prevalence = prevalence,
      as_proportion = as_proportion
    )
  )

}

#' @export
print.epiworld_virus <- function(x, ...) {
  invisible(print_virus_cpp(x))
}


#' @export
#' @details
#' The name of the `epiworld_virus` object can be manipulated with the functions
#' [set_name_virus()] and [get_name_virus()].
#' @returns
#' - The `set_name_virus` function does not return a value, but merely assigns
#' a name to the virus of choice.
#' @rdname virus
set_name_virus <- function(virus, name) {
  stopifnot_virus(virus)
  invisible(set_name_virus_cpp(virus, name))
}

#' @export
#' @returns
#' - The `get_name_virus` function returns the name of the virus of class
#' [epiworld_virus].
#' @rdname virus
get_name_virus <- function(virus) {
  stopifnot_virus(virus)
  get_name_virus_cpp(virus)
}

# Virus add --------------------------------------------------------------------

#' @export
#' @rdname virus
#' @param model An object of class `epiworld_model`.
#' @param virus An object of class `epiworld_virus`
#' @param proportion Deprecated.
#' @returns
#' - The `add_virus` function does not return a value, instead it adds the
#' virus of choice to the model object of class [epiworld_model].
add_virus <- function(model, virus, proportion) {

  if (!missing(proportion)) {

    warning(
      "The argument 'proportion' is deprecated and will be removed in ",
      "the next version."
    )

    set_distribution_virus(
      virus = virus,
      distfun = distribute_virus_randomly(proportion, as_proportion = TRUE)
    )

  } else if (isTRUE(attr(tool, "uses_deprecated"))) {

    set_distribution_virus(
      virus = tool,
      distfun = distribute_virus_randomly(
        prevalence = attr(tool, "deprecated_args")$prevalence,
        as_proportion = attr(tool, "deprecated_args")$as_proportion
      )
    )

  }

  UseMethod("add_virus")

}

#' @export
add_virus.epiworld_model <- function(model, virus, proportion) {

  stopifnot_virus(virus)
  add_virus_cpp(model, virus)
  invisible(model)

}

#' @export
add_virus.epiworld_sir <- function(model, virus, proportion) {

  stopifnot_virus(virus)
  virus_set_state(virus, init = 1, end = 2, removed = 2)
  invisible(add_virus_cpp(model, virus))

}

#' @export
add_virus.epiworld_sird <- function(model, virus, proportion) {

  stopifnot_virus(virus)
  virus_set_state(virus, init = 1, end = 2, removed = 3)
  invisible(add_virus_cpp(model, virus))

}

#' @export
add_virus.epiworld_sirconn <- function(model, virus, proportion) {

  stopifnot_virus(virus)
  add_virus.epiworld_sir(model, virus)

}

#' @export
add_virus.epiworld_sirdconn <- function(model, virus, proportion) {

  stopifnot_virus(virus)
  add_virus.epiworld_sird(model, virus)

}


#' @export
add_virus.epiworld_seir <- function(model, virus, proportion) {

  stopifnot_virus(virus)
  virus_set_state(virus, init = 1, end = 3, removed = 3)
  invisible(add_virus_cpp(model, virus))

}

#' @export
add_virus.epiworld_seird <- function(model, virus, proportion) {

  stopifnot_virus(virus)
  virus_set_state(virus, init = 1, end = 3, removed = 4)
  invisible(add_virus_cpp(model, virus))

}

#' @export
add_virus.epiworld_seirconn <- function(model, virus, proportion) {

  stopifnot_virus(virus)
  add_virus.epiworld_seir(model, virus)

}

#' @export
add_virus.epiworld_seirdconn <- function(model, virus, proportion) {

  stopifnot_virus(virus)
  add_virus.epiworld_seird(model, virus)

}

# Virus MISC -------------------------------------------------------------------

#' @export
#' @rdname virus
#' @param init,end,removed states after acquiring a virus, removing a virus,
#' and removing the agent as a result of the virus, respectively.
#' @returns
#' - The `virus_set_state` function does not return a value but assigns
#' epidemiological properties to the specified virus of class [epiworld_virus].
virus_set_state <- function(virus, init, end, removed) {

  stopifnot_virus(virus)
  invisible(virus_set_state_cpp(virus, init, end, removed))

}

#' @export
#' @returns
#' - The `rm_virus` function does not return a value, but instead removes
#' a specified virus from the model of class [epiworld_model].
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
#'   transmission_rate = 0.9, recovery = 0.1
#' )
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
#' hist(plogis(dat %*% rbind(-1, 1)))
#'
#' vfun # printing
#'
#' set_prob_infecting_fun(
#'   virus = get_virus(sir, 0),
#'   model = sir,
#'   vfun  = vfun
#' )
#'
#' run(sir, ndays = 50, seed = 11)
#' plot(sir)
#'
virus_fun_logit <- function(vars, coefs, model) {

  stopifnot_model(model)

  structure(
    virus_fun_logit_cpp(vars, coefs, model),
    class   = "epiworld_virus_fun",
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
#' @returns
#' - The `set_prob_infecting` function does not return a value, but instead
#' assigns a probability to infection for the specified virus of class
#' [epiworld_virus].
#' @rdname virus
set_prob_infecting <- function(virus, prob) {

  stopifnot_virus(virus)
  invisible(set_prob_infecting_cpp(virus, as.numeric(prob)))

}

#' @export
#' @param param Character scalar. Name of the parameter featured in `model` that
#' will be added to the virus (see details).
#' @details
#' In the case of `set_prob_infecting_ptr`, `set_prob_recovery_ptr`, and
#' `set_prob_death_ptr`, the corresponding parameters is passed as a pointer to
#' the virus. The implication of using pointers is that the values will be
#' read directly from the `model` object, so changes will be reflected.
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
#' @returns
#' - The `set_prob_recovery` function does not return a value, but instead
#' assigns a probability to recovery for the specified virus of class
#' [epiworld_virus].
#' @rdname virus
set_prob_recovery <- function(virus, prob) {

  stopifnot_virus(virus)
  invisible(set_prob_recovery_cpp(virus, as.numeric(prob)))

}

#' @export
#' @rdname virus
set_prob_recovery_ptr <- function(virus, model, param) {

  stopifnot_virus(virus)
  stopifnot_model(model)
  invisible(set_prob_recovery_ptr_cpp(virus, model, param))

}

#' @export
#' @rdname virus
set_prob_recovery_fun <- function(virus, model, vfun) {

  stopifnot_virus(virus)
  stopifnot_model(model)
  stopifnot_vfun(vfun)
  invisible(set_prob_recovery_fun_cpp(virus, model, vfun))

}

#' @export
#' @returns
#' - The `set_prob_death` function does not return a value, but instead
#' assigns a probability to death for the specified virus of class
#' [epiworld_virus].
#' @rdname virus
set_prob_death <- function(virus, prob) {

  stopifnot_virus(virus)
  invisible(set_prob_death_cpp(virus, as.numeric(prob)))

}

#' @export
#' @rdname virus
set_prob_death_ptr <- function(virus, model, param) {

  stopifnot_virus(virus)
  stopifnot_model(model)
  invisible(set_prob_death_ptr_cpp(virus, model, param))

}

#' @export
#' @rdname virus
set_prob_death_fun <- function(virus, model, vfun) {

  stopifnot_virus(virus)
  stopifnot_model(model)
  stopifnot_vfun(vfun)
  invisible(set_prob_death_fun_cpp(virus, model, vfun))

}

#' @export
#' @return
#' - The `set_incubation` function does not return a value, but instead
#' assigns an incubation period to the specified virus of class [epiworld_virus].
#' @rdname virus
set_incubation <- function(virus, incubation) {

  stopifnot_virus(virus)
  invisible(set_incubation_cpp(virus, as.numeric(incubation)))

}

#' @export
#' @rdname virus
set_incubation_ptr <- function(virus, model, param) {

  stopifnot_virus(virus)
  stopifnot_model(model)
  invisible(set_incubation_ptr_cpp(virus, model, param))

}

#' @export
#' @rdname virus
set_incubation_fun <- function(virus, model, vfun) {

  stopifnot_virus(virus)
  stopifnot_model(model)
  stopifnot_vfun(vfun)
  invisible(set_incubation_fun_cpp(virus, model, vfun))

}

#' @export
#' @rdname virus
#' @param distfun An object of class `epiworld_distribution_virus`.
set_distribution_virus <- function(virus, distfun) {

  stopifnot_virus(virus)
  stopifnot_virus_distfun(distfun)
  invisible(set_distribution_virus_cpp(virus, distfun))

}

#' @export
#' @rdname virus
#' @details The `distribute_virus_randomly` function is a factory function
#' used to randomly distribute the virus in the model. The prevalence can be set
#' as a proportion or as a number of agents. The resulting function can then be
#' passed to `set_distribution_virus`.
#' @param prevalence Numeric scalar. Prevalence of the virus.
#' @param as_proportion Logical scalar. If `TRUE`, the prevalence is set as a
#' proportion of the total number of agents in the model.
#' @return
#' - The `distribute_virus_randomly` function returns a function that can be
#' used to distribute the virus in the model. When `agents_ids` is not empty,
#' it will distribute the virus randomly within that set. Otherwise it uses
#' all the agents in the model.
distribute_virus_randomly <- function(
  prevalence,
  as_proportion,
  agents_ids = integer(0)
) {

  structure(
    distribute_virus_randomly_cpp(
      as.double(prevalence),
      as.logical(as_proportion),
      as.integer(agents_ids)
    ),
    class = "epiworld_virus_distfun"
  )

}

#' @export
#' @rdname virus
#' @param agents_ids Integer vector. Indices of the agents that will receive the
#' virus.
distribute_virus_to_set <- function(agents_ids) {

  structure(
    distribute_virus_to_set_cpp(as.vector(agents_ids)),
    class = "epiworld_virus_distfun"
  )

}

#' @export
#' @rdname virus
distribute_virus_set <- function(agents_ids) {

  .Deprecated("distribute_virus_to_set")

}
