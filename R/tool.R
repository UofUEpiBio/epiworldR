#' Tools in epiworld
#' @param model Model
#' @param name Name of the tool
#' @param susceptibility_reduction Numeric. Proportion it reduces susceptibility.
#' @param transmission_reduction Numeric. Proportion it reduces transmission.
#' @param recovery_enhancer Numeric. Proportion it improves recovery.
#' @param death_reduction Numeric. Proportion it reduces probability of death.e
#' @param tool_pos Positive integer. Index of the tool's position in the model.
#' @examples 
#' # Simple model
#' model_sirconn <- ModelSIRCONN(
#'   name                = "COVID-19",
#'   n                   = 10000,
#'   prevalence          = 0.01,
#'   contact_rate        = 5,
#'   prob_transmission   = 0.4,
#'   prob_recovery       = 0.95
#' )
#' 
#' # Running and printing
#' run(model_sirconn, ndays = 100, seed = 1912)
#' plot(model_sirconn)
#' 
#' epitool <- tool(
#'   name = "Vaccine",
#'   susceptibility_reduction = .9,
#'   transmission_reduction = .5,
#'   recovery_enhancer = .5, 
#'   death_reduction = .9
#' )
#' 
#' epitool
#' 
#' add_tool(model_sirconn, epitool, .5)
#' run(model_sirconn, ndays = 100, seed = 1912)
#' model_sirconn
#' plot(model_sirconn)
#' 
#' 
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

#' @export
print.epiworld_tool <- function(x, ...) {
  invisible(print_tool_cpp(x))
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
#' @details
#' The name of the `epiworld_tool` object can be manipulated with the functions
#' [set_name_tool()] and [get_name_tool()].
#' 
#' @rdname tool
set_name_tool <- function(tool, name) {
  stopifnot_tool(tool)
  invisible(set_name_tool_cpp(tool, name))
}

#' @export
#' @rdname tool
get_name_tool <- function(tool) {
  stopifnot_tool(tool)
  get_name_tool(tool)
}

#' @export
#' @param tool An object of class `epiworld_tool`
#' @param prevalence In the case of `add_tool`, a proportion, otherwise, an integer.
#' @rdname tool
add_tool <- function(model, tool, prevalence) UseMethod("add_tool")

#' @export
add_tool.epiworld_model <- function(model, tool, prevalence) {
  add_tool_cpp(model, tool, prevalence)
  invisible(model)
}

#' @export
#' @rdname tool
add_tool_n <- function(model, tool, prevalence) UseMethod("add_tool_n")

#' @export
add_tool_n.epiworld_model <- function(model, tool, prevalence) {
  add_tool_n_cpp(model, tool, prevalence)
  invisible(model)
}

#' @export
#' @rdname tool
rm_tool <- function(model, tool_pos) {
  invisible(rm_tool_cpp(model, tool_pos))
}

#' @export
#' @rdname tool
rm_tool <- function(model, tool_pos) {
  invisible(rm_tool_cpp(model, tool_pos))
}

# Tool functions ---------------------------------------------------------------

#' @export
#' @param vars Integer vector. Indices (starting from 0) of the positions of the
#' variables used to compute the logit probability.
#' @param coefs Numeric vector. Of the same length of `vars`, is a vector of
#' coefficients associated to the logit probability.
#' @rdname tool
#' @examples
#' 
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
#' # Creating a tool
#' mask_wearing <- tool(
#'   name = "Mask",
#'   susceptibility_reduction = 0.0,
#'   transmission_reduction   = 0.3, # Only transmission
#'   recovery_enhancer        = 0.0,
#'   death_reduction          = 0.0
#' )
#' 
#' add_tool(sir, mask_wearing, .5)
#' 
#' run(sir, ndays = 50, seed = 11)
#' hist_0 <- get_hist_total(sir)
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
#' tfun <- tool_fun_logit(
#'   vars  = c(0L, 1L),
#'   coefs = c(-1, 1),
#'   model = sir
#' )
#' 
#' # The infection prob is lower
#' hist(plogis(dat %*% rbind(.5,1)))
#' 
#' tfun # printing
#' 
#' 
#' set_susceptibility_reduction_fun(
#'   tool  = get_tool(sir, 0),
#'   model = sir,
#'   tfun  = tfun
#'   )
#'   
#' run(sir, ndays = 50, seed = 11)
#' hist_1 <- get_hist_total(sir)
#' 
#' op <- par(mfrow = c(1, 2))
#' plot(hist_0); abline(v = 30)
#' plot(hist_1); abline(v = 30)
#' par(op)
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

# Susceptibility reduction -----------------------------------------------------


#' @export
#' @param prob Numeric scalar. A probability (between zero and one).
#' @rdname tool
set_susceptibility_reduction <- function(tool, prob) {
  
  stopifnot_tool(tool)
  set_susceptibility_reduction_cpp(tool, prob)
  
}

#' @export
#' @param param Character scalar. Name of the parameter featured in `model` that
#' will be added to the tool (see details.)
#' @details
#' In the case of `set_susceptibility_reduction_ptr`, `set_transmission_reduction_ptr`, 
#' `set_recovery_enhancer`, and
#' `set_death_reduction_ptr`, the corresponding parameters is passed as a pointer to
#' the tool. The implication of using pointers is that the values will be
#' read directly from the `model` object, so changes will be reflected.
#' 
#' @rdname tool
set_susceptibility_reduction_ptr <- function(tool, model, param) {
  
  stopifnot_tool(tool)
  stopifnot_model(model)
  invisible(set_susceptibility_reduction_ptr_cpp(tool, model, param))
  
}

#' @export
#' @param tfun An object of class `epiworld_tool_fun`.
#' @rdname tool
set_susceptibility_reduction_fun <- function(tool, model, tfun) {
  
  stopifnot_tool(tool)
  stopifnot_model(model)
  stopifnot_tfun(tfun)
  invisible(set_susceptibility_reduction_fun_cpp(tool, model, tfun))
  
}

# Transmission reduction -------------------------------------------------------

#' @export
#' @rdname tool
set_transmission_reduction <- function(tool, prob) {
  
  stopifnot_tool(tool)
  set_transmission_reduction_cpp(tool, prob)
  
}

#' @export
#' @rdname tool
set_transmission_reduction_ptr <- function(tool, model, param) {
  
  stopifnot_tool(tool)
  stopifnot_model(model)
  set_transmission_reduction_ptr_cpp(tool, model, param)
  
}

#' @export
#' @rdname tool
set_transmission_reduction_fun <- function(tool, model, tfun) {
  
  stopifnot_tool(tool)
  stopifnot_model(model)
  stopifnot_tfun(tfun)
  set_transmission_reduction_fun_cpp(tool, model, tfun)
}

# Recovery enhancer ------------------------------------------------------------

#' @export
#' @rdname tool
set_recovery_enhancer <- function(tool, prob) {
  
  stopifnot_tool(tool)
  set_recovery_enhancer_cpp(tool, prob)
  
}

#' @export
#' @rdname tool
set_recovery_enhancer_ptr <- function(tool, model, param) {
  
  stopifnot_tool(tool)
  stopifnot_model(model)
  set_recovery_enhancer_ptr_cpp(tool, model, param)
  
}

#' @export
#' @rdname tool
set_recovery_enhancer_fun <- function(tool, model, tfun) {
  
  stopifnot_tool(tool)
  stopifnot_model(model)
  stopifnot_tfun(tfun)
  set_recovery_enhancer_fun_cpp(tool, model, tfun)
  
}

# Death reduction --------------------------------------------------------------

#' @export
#' @rdname tool
set_death_reduction <- function(tool, prob) {
  
  stopifnot_tool(tool)
  set_death_reduction_cpp(tool, prob)
  
}

#' @export
#' @rdname tool
set_death_reduction_ptr <- function(tool, model, param) {
  
  stopifnot_tool(tool)
  stopifnot_model(model)
  set_death_reduction_ptr_cpp(tool, model, param)
  
}

#' @export
#' @rdname tool
set_death_reduction_fun <- function(tool, model, tfun) {
  
  stopifnot_tool(tool)
  stopifnot_model(model)
  stopifnot_tfun(tfun)
  set_death_reduction_fun_cpp(tool, model, tfun)
  
}


