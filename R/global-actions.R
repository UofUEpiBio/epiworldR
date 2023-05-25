
#' Global Actions
#' 
#' Global actions are functions that are executed at each time step of the
#' simulation. They are useful for implementing interventions, such as
#' vaccination, isolation, and social distancing by means of tools.
#' 
#' @export
#' @param prob Numeric scalar. A probability between 0 and 1.
#' @param tool An object of class [tool].
#' @name global-actions
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
#' # Creating a tool
#' epitool <- tool(
#'   name = "Vaccine",
#'   susceptibility_reduction = .9,
#'   transmission_reduction = .5,
#'   recovery_enhancer = .5, 
#'   death_reduction = .9
#' )
#' 
#' 
#' # Adding a global action
#' vaccine_day_20 <- globalaction_tool(epitool, .2, day = 20)
#' add_global_action(model_sirconn, vaccine_day_20)
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
#'   prob_transmission   = 0.4,
#'   prob_recovery       = 0.95
#' )
#' 
#' closure_day_10 <- globalaction_set_params("Contact rate", 0, day = 10)
#' add_global_action(model_sirconn2, closure_day_10)
#' 
#' # Running and printing
#' run(model_sirconn2, ndays = 40, seed = 1912)
#' model_sirconn2
#' plot_incidence(model_sirconn2)
globalaction_tool <- function(
  tool, prob,
  name = get_name_tool(tool), day = -99
  ) {

  structure(
    globalaction_tool_cpp(tool, prob, name, day),
    class = c("epiworld_globalaction_tool", "epiworld_globalaction"),
    tool = tool,
    call = match.call()
  )
}

#' @export
#' @rdname global-actions
#' @param vars Integer vector. The position of the variables in the model.
#' @param coefs Numeric vector. The coefficients of the logistic regression.
#' @details The function `globalaction_tool_logit` allows to specify a logistic
#' regression model for the probability of using a tool. The model is specified
#' by the vector of coefficients `coefs` and the vector of variables `vars`.
#' `vars` is an integer vector indicating the position of the variables in the
#' model.
globalaction_tool_logit <- function(
  tool, vars, coefs,
  name = get_name_tool(tool), day = -99
  ) {

  structure(
    globalaction_tool_logit_cpp(tool, vars, coefs, name, day),
    class = c("epiworld_globalaction_tool_logit", "epiworld_globalaction"),
    tool = tool,
    call = match.call()
  )

}

#' @export 
#' @param param Character scalar. The name of the parameter to be set.
#' @param value Numeric scalar. The value of the parameter.
#' @rdname global-actions
#' @details The function `globalaction_set_param` allows to set a parameter of
#' the model. The parameter is specified by its name `param` and the value by
#' `value`.
globalaction_set_params <- function(
  param, value,
  name = paste0("Set ", param, " to ", value), day = -99
  ) {

  structure(
    globalaction_set_param_cpp(param, value, name, day),
    class = c("epiworld_globalaction_set_param", "epiworld_globalaction"),
    param = param,
    value = value
  )
}

#' @export
print.epiworld_globalaction <- function(x, ...) {

  print_global_action_cpp(x)
  cat("Call: ", deparse(attr(x, "call")), "\n")
  if (length(attr(x, "tool"))) {
    cat("Tool: ", get_name_tool(attr(x, "tool")), "\n")
  } else if (inherits(x, "epiworld_globalaction_set_param")) {
    cat("Parameter: ", attr(x, "param"), "\n")
    cat("Value: ", attr(x, "value"), "\n")
  }

  invisible(x)

}

#' @export
#' @param action A global action.
#' @param date Integer. The date at which the action is executed (see details).
#' @param model An object of class [epiworld_model].
#' @param name Character scalar. The name of the action.
#' @rdname global-actions
#' @seealso epiworld-model
#' @details The function `add_global_action` adds a global action to a model.
#' The model checks for actions to be executed at each time step. If the added
#' action matches the current time step, the action is executed. When `date` is
#' negative, the action is executed at each time step. When `date` is positive,
#' the action is executed at the specified time step.
add_global_action <- function(model, action) {
  
  if (length(attr(action, "tool")))
    add_tool_n(model, attr(action, "tool"), 0)

  invisible(add_global_action_cpp(model, action))

}

