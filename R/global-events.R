
#' Global Events
#'
#' Global events are functions that are executed at each time step of the
#' simulation. They are useful for implementing interventions, such as
#' vaccination, isolation, and social distancing by means of tools.
#'
#' @export
#' @param prob Numeric scalar. A probability between 0 and 1.
#' @param tool An object of class [tool].
#' @param day Integer. The day (step) at which the event is executed (see details).
#' When negative (default -99), the event is executed at each time step.
#' When positive, the event is executed only at the specified day.
#' @param name Character scalar. The name of the event.
#' @name global-events
#' @aliases global-actions
#' @details 
#' Global events can be created using:
#' - `globalevent_tool()`: Distribute a tool to agents with a given probability
#' - `globalevent_tool_logit()`: Distribute a tool using a logistic regression model
#' - `globalevent_set_params()`: Change a model parameter at a specific time
#' - `globalevent_fun()`: Execute a custom R function at each time step or specific day
#' 
#' Once created, events are added to a model with `add_globalevent()` and can be 
#' removed by name with `rm_globalevent()`.
#' 
#' **Note:** The `globalaction_*` family of functions are deprecated. Use the 
#' `globalevent_*` functions instead.
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
#' # Creating a tool
#' epitool <- tool(
#'   name = "Vaccine",
#'   prevalence = 0,
#'   as_proportion = FALSE,
#'   susceptibility_reduction = .9,
#'   transmission_reduction = .5,
#'   recovery_enhancer = .5,
#'   death_reduction = .9
#' )
#'
#'
#' # Adding a global event
#' vaccine_day_20 <- globalevent_tool(epitool, .2, day = 20)
#' add_globalevent(model_sirconn, vaccine_day_20)
#'
#' # Running and printing
#' run(model_sirconn, ndays = 40, seed = 1912)
#' model_sirconn
#' plot_incidence(model_sirconn)
#'
#' # Example 2: Changing the contact rate -------------------------------------
#' model_sirconn2 <- ModelSIRCONN(
#'   name                = "COVID-19",
#'   n                   = 10000,
#'   prevalence          = 0.01,
#'   contact_rate        = 5,
#'   transmission_rate   = 0.4,
#'   recovery_rate       = 0.95
#' )
#'
#' closure_day_10 <- globalevent_set_params("Contact rate", 0, day = 10)
#' add_globalevent(model_sirconn2, closure_day_10)
#'
#' # Running and printing
#' run(model_sirconn2, ndays = 40, seed = 1912)
#' model_sirconn2
#' plot_incidence(model_sirconn2)
#' @returns
#' - The `globalevent_set_params` function returns an object of class
#' [epiworld_globalevent_set_param] and [epiworld_globalevent].
#'
#' - `globalevent_tool` returns an object of class
#' [epiworld_globalevent_tool] and [epiworld_globalevent].
#'
#' - `globalevent_tool_logit` returns an object of class
#' [epiworld_globalevent_tool_logit] and [epiworld_globalevent].
#' @aliases
#' epiworld_globalevent_set_param
#' epiworld_globalevent_tool
#' epiworld_globalevent_tool_logit
#' epiworld_globalevent
#' actions
#'
globalevent_tool <- function(
    tool, prob,
    name = get_name_tool(tool), day = -99
    ) {

  structure(
    globalevent_tool_cpp(tool, prob, name, day),
    class = c("epiworld_globalevent_tool", "epiworld_globalevent"),
    tool = tool,
    call = match.call()
  )
}

#' @export
#' @rdname epiworldR-deprecated
globalaction_tool <- function(...) {

  .Defunct(
    new = "globalevent_tool"
  )

}

#' @export
#' @rdname global-events
#' @param vars Integer vector. The position of the variables in the model.
#' @param coefs Numeric vector. The coefficients of the logistic regression.
#' @details The function `globalevent_tool_logit` allows to specify a logistic
#' regression model for the probability of using a tool. The model is specified
#' by the vector of coefficients `coefs` and the vector of variables `vars`.
#' `vars` is an integer vector indicating the position of the variables in the
#' model.
globalevent_tool_logit <- function(
    tool, vars, coefs,
    name = get_name_tool(tool), day = -99
    ) {

  stopifnot_tool(tool)

  structure(
    globalevent_tool_logit_cpp(
      tool,
      as.integer(vars),
      as.double(coefs),
      name,
      as.integer(day)
    ),
    class = c("epiworld_globalevent_tool_logit", "epiworld_globalevent"),
    tool = tool,
    call = match.call()
  )

}

#' @export
#' @rdname epiworldR-deprecated
globalaction_tool_logit <- function(...) {

  .Defunct(
    new = "globalevent_tool_logit"
  )

  globalevent_tool_logit(...)

}

#' @export
#' @param param Character scalar. The name of the parameter to be set.
#' @param value Numeric scalar. The value of the parameter.
#' @rdname global-events
#' @details The function `globalevent_set_param` allows to set a parameter of
#' the model. The parameter is specified by its name `param` and the value by
#' `value`.
globalevent_set_params <- function(
    param, value,
    name = paste0("Set ", param, " to ", value), day = -99
    ) {

  structure(
    globalevent_set_param_cpp(
      param,
      as.double(value),
      name,
      as.integer(day)
    ),
    class = c("epiworld_globalevent_set_param", "epiworld_globalevent"),
    param = param,
    value = as.double(value),
    call = match.call()
  )
}

#' @export
#' @rdname epiworldR-deprecated
globalaction_set_params <- function(...) {


  .Defunct(
    new = "globalevent_set_params"
  )

  globalevent_set_params(...)

}

#' @export
#' @rdname global-events
#' @param fun Function. The function to be executed.
#' @details The function `globalevent_fun` allows to specify a function to be
#' executed at a given day. The function object must receive an object of class
#' [epiworld_model] as only argument.
#' @examples
#' # Example using `globalevent_fun` to record the state of the
#' # agents at each time step.
#'
#' # We start by creating an SIR connected model
#' model <- ModelSIRCONN(
#'   name              = "SIR with Global Saver",
#'   n                 = 1000,
#'   prevalence        = 0.01,
#'   contact_rate      = 5,
#'   transmission_rate = 0.4,
#'   recovery_rate     = 0.3
#' )
#'
#' # We create the object where the history of the agents will be stored
#' agents_history <- NULL
#'
#' # This function prints the total number of agents in each state
#' # and stores the history of the agents in the object `agents_history`
#' hist_saver <- function(m) {
#'
#'   message("Today's totals are: ", paste(get_today_total(m), collapse = ", "))
#'
#'   # We use the `<<-` operator to assign the value to the global variable
#'   # `agents_history` (see ?"<<-")
#'   agents_history <<- cbind(
#'     agents_history,
#'     get_agents_states(m)
#'   )
#'
#' }
#'
#' # We create the global event that will execute the function `hist_saver`
#' # at each time step
#' hist_saver_event <- globalevent_fun(hist_saver, "Agent History Saver")
#'
#' # We add the global event to the model
#' model <- add_globalevent(model, hist_saver_event)
globalevent_fun <- function(
    fun, name = deparse(substitute(fun)), day = -99
    ) {

  structure(
    globalevent_fun_cpp(fun, name, as.integer(day)),
    class = c("epiworld_globalevent_fun", "epiworld_globalevent"),
    fun = fun,
    call = match.call()
  )

}

#' @export
#' @rdname epiworldR-deprecated
globalaction_fun <- function(...) {

  .Defunct(
    new = "globalevent_fun"
  )

  globalevent_fun(...)

}

#' @export
print.epiworld_globalevent <- function(x, ...) {

  print_global_action_cpp(x)
  cat("Call: ", deparse(attr(x, "call")), "\n")
  if (length(attr(x, "tool"))) {
    cat("Tool: ", get_name_tool(attr(x, "tool")), "\n")
  } else if (inherits(x, "epiworld_globalevent_set_param")) {
    cat("Parameter: ", attr(x, "param"), "\n")
    cat("Value: ", attr(x, "value"), "\n")
  }

  invisible(x)

}

#' Add Global Event to Model
#'
#' @export
#' @param model An object of class [epiworld_model].
#' @param event An object of class `epiworld_globalevent` to be added to the model.
#' @param action (Deprecated) use `event` instead.
#' @rdname add_globalevent
#' @seealso epiworld-model
#' @details The function `add_globalevent` adds a global event to a model.
#' The model checks for events to be executed at each time step. If the added
#' event matches the current time step, the event is executed. When `day` is
#' negative (the default), the event is executed at each time step. When `day` 
#' is positive, the event is executed at the specified time step.
#' @returns The model with the added event (invisibly).
#' @examples
#' # See examples in ?global-events
add_globalevent <- function(model, event, action = NULL) {

  if (missing(event) && !missing(action)) {
    event <- action
    warning("The argument `action` is deprecated. Use `event` instead.")
  }

  stopifnot_model(model)

  if (length(attr(event, "tool")))
    add_tool(model, attr(event, "tool"))

  invisible(add_globalevent_cpp(model, event))

}

#' Remove Global Event from Model
#'
#' @export
#' @param model An object of class [epiworld_model].
#' @param name Character scalar. The name of the global event to remove.
#' @rdname rm_globalevent
#' @seealso epiworld-model
#' @details The function `rm_globalevent` removes a global event from a model
#' by its name. The name should match the name assigned when the event was 
#' created (e.g., through `globalevent_tool()`, `globalevent_fun()`, etc.).
#' @returns The model with the event removed (invisibly).
#' @examples
#' \dontrun{
#' # Create a model
#' model_sirconn <- ModelSIRCONN(
#'   name = "COVID-19",
#'   n = 10000,
#'   prevalence = 0.01,
#'   contact_rate = 5,
#'   transmission_rate = 0.4,
#'   recovery_rate = 0.95
#' )
#'
#' # Add a global event
#' epitool <- tool(
#'   name = "Vaccine",
#'   prevalence = 0,
#'   as_proportion = FALSE,
#'   susceptibility_reduction = .9,
#'   transmission_reduction = .5,
#'   recovery_enhancer = .5,
#'   death_reduction = .9
#' )
#' vaccine_event <- globalevent_tool(epitool, .2, name = "Vaccination", day = 20)
#' add_globalevent(model_sirconn, vaccine_event)
#'
#' # Remove the global event by name
#' rm_globalevent(model_sirconn, "Vaccination")
#' }
rm_globalevent <- function(model, name) {

  stopifnot_model(model)
  invisible(rm_globalevent_cpp(model, name))

}
