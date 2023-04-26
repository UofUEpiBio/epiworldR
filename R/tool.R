#' Tools in epiworld
#' @param m Model
#' @param name Name of the tool
#' @param susceptibility_reduction Numeric. Proportion it reduces susceptibility.
#' @param transmission_reduction Numeric. Proportion it reduces transmission.
#' @param recovery_enhancer Numeric. Proportion it improves recovery.
#' @param death_reduction Numeric. Proportion it reduces probability of death.e
#' @param tool_pos Positive integer. Index of the tool's position in the model.
#' @examples 
#' epitool <- tool(
#'   name = "Vaccine",
#'   susceptibility_reduction = .9,
#'   transmission_reduction = .5,
#'   recovery_enhancer = .5, 
#'   death_reduction = .9
#' )
#' @export
tool <- function(
    name,
    susceptibility_reduction,
    transmission_reduction,
    recovery_enhancer,
    death_reduction
) {

  structure(
    tool_cpp(
      name,
      susceptibility_reduction,
      transmission_reduction,
      recovery_enhancer,
      death_reduction
    ),
    class = "epiworld_tool"
  )
    
}

stopifnot_tool <- function(tool) {
  if (!inherits(tool, "epiworld_tool")) {
    stop(
      "The -tool- object must be of class \"epiworld_tool\". ",
      "The object passed to the function is of class(es): ", 
      paste(class(tool), collapse = ", ")
    )
  }
}

stopifnot_tfun <- function(tfun) {
  if (!inherits(tfun, "epiworld_tool_fun")) {
    stop(
      "The -tfun- object must be of class \"epiworld_tool_fun\". ",
      "The object passed to the function is of class(es): ", 
      paste(class(tfun), collapse = ", ")
    )
  }
}

#' @export
#' @param t An object of class `epiworld_tool`
#' @param prevalence In the case of `add_tool`, a proportion, otherwise, an integer.
#' @rdname tool
add_tool <- function(m, t, prevalence) UseMethod("add_tool")

#' @export
add_tool.epiworld_model <- function(m, t, prevalence) {
  add_tool_cpp(m, t, prevalence)
  invisible(m)
}

#' @export
#' @rdname tool
add_tool_n <- function(m, t, prevalence) UseMethod("add_tool_n")

#' @export
add_tool_n.epiworld_model <- function(m, t, prevalence) {
  add_tool_n_cpp(m, t, prevalence)
  invisible(m)
}

#' @export
#' @rdname tool
rm_tool <- function(m, tool_pos) {
  invisible(rm_tool_cpp(m, tool_pos))
}

#' @export
#' @rdname tool
rm_tool <- function(m, tool_pos) {
  invisible(rm_tool_cpp(m, tool_pos))
}

# Tool functions ---------------------------------------------------------------

#' @export
#' @param vars Integer vector. Indices (starting from 0) of the positions of the
#' variables used to compute the logit probability.
#' @param coefs Numeric vector. Of the same length of `vars`, is a vector of
#' coefficients associated to the logit probability.
#' @rdname tool
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
#' set_susceptibility_reductionfun(
#'   virus = get_virus(sir, 0),
#'   model = sir,
#'   vfun  = vfun
#'   )
#'   
#' run(sir, ndays = 50, seed = 11)
#' plot(sir)
#' 
#' 
tool_fun_logit <- function(vars, coefs, model) {
  
  stopifnot_model(model)
  
  structure(
    tool_fun_logit_cpp(vars, coefs, model),
    class = "epiworld_tool_fun",
    builder = "tool_fun_logit",
    vars    = vars,
    coefs   = coefs,
    model   = model
  )
  
}

#' @export
print.epiworld_tool_fun <- function(x, ...) {
  
  cat("An epiworld_tool_function object.\n")
  cat("(model: ", get_name(attr(x, "model")), ")\n", sep = "")
  cat("This function was built using -tool_fun_logit()-. and it features ")
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
#' @rdname tool
set_susceptibility_reduction <- function(tool, prob) {
  
  stopifnot_virus(tool)
  set_susceptibility_reduction_cpp(tool, prob)
  
}

#' @export
#' @param param Character scalar. Name of the parameter featured in `model` that
#' will be added to the virus (see details.)
#' @details
#' In the case of `set_susceptibility_reductionptr`, `set_transmission_reduction_ptr`, and
#' `set_prob_death_ptr`, the corresponding parameters is passed as a pointer to
#' the virus. The implication of using pointers is that the values will be
#' read directly from the `model` object, so changes will be reflected.
#' 
#' @rdname tool
set_susceptibility_reduction_ptr <- function(tool, model, param) {
  
  stopifnot_tool(tool)
  stopifnot_model(model)
  invisible(set_susceptibility_reduction_ptr_cpp(tool, model, param))
  
}

#' @export
#' @rdname tool
set_susceptibility_reduction_fun <- function(tool, model, tfun) {
  
  stopifnot_tool(tool)
  stopifnot_model(model)
  stopifnot_tfun(tfun)
  invisible(set_susceptibility_reductionfun_cpp(tool, model, tfun))
  
}

#' @export
#' @rdname tool
set_transmission_reduction <- function(virus, prob) {
  
  stopifnot_tool(tool)
  set_transmission_reduction_cpp(virus, prob)
  
}

#' @export
#' @rdname virus
set_transmission_reduction_ptr <- function(virus, model, param) {
  
  stopifnot_tool(tool)
  stopifnot_model(model)
  set_transmission_reduction_ptr_cpp(virus, model, param)
  
}

#' @export
#' @rdname virus
set_transmission_reduction_fun <- function(virus, model, vfun) {
  
  stopifnot_tool(tool)
  stopifnot_model(model)
  stopifnot_tfun(vfun)
  set_transmission_reduction_fun_cpp(virus, model, vfun)
}

#' @export
#' @rdname virus
set_prob_death <- function(virus, prob) {
  
  stopifnot_tool(tool)
  set_prob_death_cpp(virus, prob)
  
}

#' @export
#' @rdname virus
set_prob_death_ptr <- function(virus, model, param) {
  
  stopifnot_tool(tool)
  stopifnot_model(model)
  set_prob_death_ptr_cpp(virus, model, param)
  
}

#' @export
#' @rdname virus
set_prob_death_fun <- function(virus, model, vfun) {
  
  stopifnot_tool(tool)
  stopifnot_model(model)
  stopifnot_tfun(vfun)
  set_prob_death_fun_cpp(virus, model, vfun)
  
}

