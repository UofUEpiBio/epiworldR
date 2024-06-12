#' Tools in epiworld
#' 
#' Tools are functions that affect how agents react to the virus. They can be
#' used to simulate the effects of vaccination, isolation, and social
#' distancing.
#' 
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
#'   transmission_rate   = 0.4,
#'   recovery_rate       = 0.95
#' )
#' 
#' # Running and printing
#' run(model_sirconn, ndays = 100, seed = 1912)
#' plot(model_sirconn)
#' 
#' epitool <- tool(
#'   name = "Vaccine",
#'   prevalence = 0.5,
#'   as_proportion = TRUE,
#'   susceptibility_reduction = .9,
#'   transmission_reduction = .5,
#'   recovery_enhancer = .5, 
#'   death_reduction = .9
#' )
#' 
#' epitool
#' 
#' set_name_tool(epitool, 'Pfizer') # Assigning name to the tool
#' get_name_tool(epitool) # Returning the name of the tool
#' add_tool(model_sirconn, epitool)
#' run(model_sirconn, ndays = 100, seed = 1912)
#' model_sirconn
#' plot(model_sirconn)
#' 
#' # To declare a certain number of individuals with the tool
#' rm_tool(model_sirconn, 0) # Removing epitool from the model
#' set_prevalence_tool(epitool, 0.1, TRUE) # Setting prevalence to 0.1
#' add_tool(model_sirconn, epitool)
#' run(model_sirconn, ndays = 100, seed = 1912)
#' 
#' # Adjusting probabilities due to tool
#' set_susceptibility_reduction(epitool, 0.1) # Susceptibility reduction 
#' set_transmission_reduction(epitool, 0.2) # Transmission reduction
#' set_recovery_enhancer(epitool, 0.15) # Probability increase of recovery
#' set_death_reduction(epitool, 0.05) # Probability reduction of death
#' 
#' rm_tool(model_sirconn, 0) 
#' add_tool(model_sirconn, epitool)
#' run(model_sirconn, ndays = 100, seed = 1912) # Run model to view changes
#' 
#' @export
#' @returns 
#' - The `tool` function creates a tool of class [epiworld_tool].
#' @aliases epiworld_tool
tool <- function(
    name,
    prevalence,
    as_proportion,
    susceptibility_reduction,
    transmission_reduction,
    recovery_enhancer,
    death_reduction
) {

  structure(
    tool_cpp(
      name,
      prevalence,
      as_proportion,
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

stopifnot_tool_distfun <- function(tool_distfun) {
  if (!inherits(tool_distfun, "epiworld_tool_distfun")) {
    stop(
      "The -tool_distfun- object must be of class \"epiworld_tool_distfun\". ",
      "The object passed to the function is of class(es): ", 
      paste(class(tool_distfun), collapse = ", ")
    )
  }
}

#' @export
#' @details
#' The name of the `epiworld_tool` object can be manipulated with the functions
#' [set_name_tool()] and [get_name_tool()].
#' @returns 
#' - The `set_name_tool` function assigns a name to the tool of class 
#' [epiworld_tool] and returns the tool.
#' @rdname tool
set_name_tool <- function(tool, name) {
  stopifnot_tool(tool)
  invisible(set_name_tool_cpp(tool, name))
}


#' @returns
#' - The `get_name_tool` function returns the name of the tool of class 
#' [epiworld_tool].
#' @rdname tool
#' @export
get_name_tool <- function(tool) {
  stopifnot_tool(tool)
  get_name_tool_cpp(tool)
}

#' @export
#' @param tool An object of class `epiworld_tool`
#' @param proportion Deprecated. Either set the prevalence during the tool
#' initialization or use `set_prevalence_tool`.
#' @details 
#' The `add_tool` function adds the specified tool to the model of class 
#' [epiworld_model] with specified proportion.
#' @rdname tool
add_tool <- function(model, tool, proportion) {
  
  if (!missing(proportion)) {

    warning(
      "The 'proportion' argument is deprecated. ",
      "Use 'set_prevalence_tool' instead."
      )

    set_prevalence_tool(tool, proportion, TRUE)

  }

  UseMethod("add_tool")

}

#' @export
add_tool.epiworld_model <- function(model, tool, proportion) {

  stopifnot_tool(tool)
  add_tool_cpp(model, tool)
  invisible(model)

}

#' @export
#' @returns 
#' - The `rm_tool` function removes the specified tool from a model.
#' @rdname tool
rm_tool <- function(model, tool_pos) {

  stopifnot_model(model)
  invisible(rm_tool_cpp(model, as.integer(tool_pos)))

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
#'   transmission_rate = 0.9, recovery_rate = 0.1
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
#'   prevalence               = 0.5,
#'   as_proportion            = TRUE,
#'   susceptibility_reduction = 0.0,
#'   transmission_reduction   = 0.3, # Only transmission
#'   recovery_enhancer        = 0.0,
#'   death_reduction          = 0.0
#' )
#' 
#' add_tool(sir, mask_wearing)
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
    tool_fun_logit_cpp(as.integer(vars), as.double(coefs), model),
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
#' @returns
#' - The `set_susceptibility_reduction` function assigns a probability reduction 
#' to the specified tool of class [epiworld_tool]. 
#' @rdname tool
set_susceptibility_reduction <- function(tool, prob) {
  
  stopifnot_tool(tool)
  set_susceptibility_reduction_cpp(tool, as.double(prob))
  
}

#' @export
#' @param param Character scalar. Name of the parameter featured in `model` that
#' will be added to the tool (see details).
#' @details
#' In the case of `set_susceptibility_reduction_ptr`, `set_transmission_reduction_ptr`, 
#' `set_recovery_enhancer`, and
#' `set_death_reduction_ptr`, the corresponding parameters are passed as a pointer to
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
#' @returns
#' - The `set_transmission_reduction` function assigns a probability reduction 
#' to the specified tool of class [epiworld_tool]. 
#' @rdname tool
set_transmission_reduction <- function(tool, prob) {
  
  stopifnot_tool(tool)
  invisible(set_transmission_reduction_cpp(tool, as.double(prob)))
  
}

#' @export
#' @rdname tool
set_transmission_reduction_ptr <- function(tool, model, param) {
  
  stopifnot_tool(tool)
  stopifnot_model(model)
  invisible(set_transmission_reduction_ptr_cpp(tool, model, param))
  
}

#' @export
#' @rdname tool
set_transmission_reduction_fun <- function(tool, model, tfun) {
  
  stopifnot_tool(tool)
  stopifnot_model(model)
  stopifnot_tfun(tfun)
  invisible(set_transmission_reduction_fun_cpp(tool, model, tfun))
}

# Recovery enhancer ------------------------------------------------------------

#' @export
#' @returns
#' - The `set_recovery_enhancer` function assigns a probability increase 
#' to the specified tool of class [epiworld_tool]. 
#' @rdname tool
set_recovery_enhancer <- function(tool, prob) {
  
  stopifnot_tool(tool)
  invisible(set_recovery_enhancer_cpp(tool, as.double(prob)))
  
}

#' @export
#' @rdname tool
set_recovery_enhancer_ptr <- function(tool, model, param) {
  
  stopifnot_tool(tool)
  stopifnot_model(model)
  invisible(set_recovery_enhancer_ptr_cpp(tool, model, param))
  
}

#' @export
#' @rdname tool
set_recovery_enhancer_fun <- function(tool, model, tfun) {
  
  stopifnot_tool(tool)
  stopifnot_model(model)
  stopifnot_tfun(tfun)
  invisible(set_recovery_enhancer_fun_cpp(tool, model, tfun))
  
}

# Death reduction --------------------------------------------------------------

#' @export
#' @returns
#' - The `set_death_reduction` function assigns a probability decrease 
#' to the specified tool of class [epiworld_tool]. 
#' @rdname tool
set_death_reduction <- function(tool, prob) {
  
  stopifnot_tool(tool)
  invisible(set_death_reduction_cpp(tool, as.double(prob)))
  
}

#' @export
#' @rdname tool
set_death_reduction_ptr <- function(tool, model, param) {
  
  stopifnot_tool(tool)
  stopifnot_model(model)
  invisible(set_death_reduction_ptr_cpp(tool, model, param))
  
}

#' @export
#' @rdname tool
set_death_reduction_fun <- function(tool, model, tfun) {
  
  stopifnot_tool(tool)
  stopifnot_model(model)
  stopifnot_tfun(tfun)
  invisible(set_death_reduction_fun_cpp(tool, model, tfun))
  
}

#' @export
#' @rdname agents_smallworld
#' @returns 
#' - `get_agents_tools` returns a list of class `epiworld_agents_tools`
#' with `epiworld_tools` (list of lists).
get_agents_tools <- function(model) {
  
  stopifnot_model(model)
  
  res <- lapply(
    get_agents_tools_cpp(model),
    `class<-`,
    "epiworld_tools"
  )
  
  structure(res, class = c("epiworld_agents_tools", class(res)))
  
}

#' @export 
#' @rdname tool
#' @param max_print Numeric scalar. Maximum number of tools to print.
#' @param ... Currently ignored.
#' @param x An object of class `epiworld_agents_tools`.
print.epiworld_agents_tools <- function(x, max_print = 10, ...) {
  
  for (i in 1:min(max_print, length(x))) {
    print_agent_tools_cpp(x[[i]])
  }
  
  if (length(x) > max_print) {
    cat(sprintf("Showing first %s of %s tools.\n", max_print, length(x)))
  }
  
  invisible(x)
  
}

#' @export 
#' @details
#' The `set_distribution_tool` function assigns a distribution function to the
#' specified tool of class [epiworld_tool]. The distribution function can be
#' created using the functions [distribute_tool_randomly()] and
#' [distribute_tool_to_set()].
#' @rdname tool
set_distribution_tool <- function(tool, distfun) {

  stopifnot_tool(tool)
  stopifnot_tool_distfun(distfun)
  invisible(set_distribution_tool_cpp(tool, distfun))

}

#' @export
#' @rdname tool
#' @details 
#' The `distribute_tool_randomly` function creates a distribution function that
#' randomly assigns the tool to a proportion of the population.
#' @return 
#' - The `distribute_tool_randomly` function returns a distribution function of
#' class `epiworld_tool_distfun`.
distribute_tool_randomly <- function(
  prevalence,
  as_proportion
) {

  structure(
    distribute_tool_randomly_cpp(
      as.double(prevalence),
      as.logical(as_proportion)
    ),
    class = "epiworld_tool_distfun"
  )

}

#' @export
#' @rdname tool
#' @details 
#' The `distribute_tool_to_set` function creates a distribution function that
#' assigns the tool to a set of agents.
#' @param agents_ids Integer vector. Indices of the agents to which the tool
#' will be assigned.
#' @return
#' - The `distribute_tool_to_set` function returns a distribution function of
#' class `epiworld_tool_distfun`.
distribute_tool_to_set <- function(
  agents_ids
) {

  structure(
    distribute_tool_to_set_cpp(
      agents_ids
    ),
    class = "epiworld_tool_distfun"
  )

}