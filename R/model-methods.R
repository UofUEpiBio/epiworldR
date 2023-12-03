stopifnot_model <- function(model) {
  if (!inherits(model, "epiworld_model")) {
    stop(
      "The -model- object must be of class \"epiworld_model\". ",
      "The object passed to the function is of class(es): ", 
      paste(class(model), collapse = ", ")
    )
  }
}

#' Methods for epiworldR objects
#' 
#' The functions described in this section are methods for objects of class
#' `epiworld_model`. Besides of printing and plotting, other methods provide
#' access to manipulate model parameters, getting information about the model
#' and running the simulation.
#' 
#' @param x An object of class `epiworld_model`.
#' @param ndays Number of days (steps) of the simulation.
#' @param seed Seed to set for initializing random number generator.
#' @param model Model object.
#' @export
#' @name epiworld-methods
#' @aliases epiworld_model
#' @examples
#' 
#' model_sirconn <- ModelSIRCONN(
#' name                = "COVID-19",
#' n                   = 10000,
#' prevalence          = 0.01,
#' contact_rate        = 5,
#' transmission_rate   = 0.4,
#' recovery_rate       = 0.95
#' )
#' 
#' # Queuing - If you wish to implement the queuing function, declare whether 
#' # you would like it "on" or "off", if any. 
#' queuing_on(model_sirconn)
#' queuing_off(model_sirconn)
#' run(model_sirconn, ndays = 100, seed = 1912)
#' 
#' # Verbose - "on" prints the progress bar on the screen while "off" 
#' # deactivates the progress bar. Declare which function you want to implement, 
#' # if any. 
#' verbose_on(model_sirconn)
#' verbose_off(model_sirconn)
#' run(model_sirconn, ndays = 100, seed = 1912)
#' 
#' get_states(model_sirconn) # Returns all unique states found within the model.
#' 
#' get_param(model_sirconn, 'Contact rate') # Returns the value of the selected 
#'                                          # parameter within the model object.
#'                                          # In order to view the parameters, 
#'                                          # run the model object and find the 
#'                                          # "Model parameters" section. 
#' 
#' set_param(model_sirconn, 'Contact rate', 2) # Allows for adjustment of model 
#'                                             # parameters within the model 
#'                                             # object. In this example, the 
#'                                             # Contact rate parameter is 
#'                                             # changed to 2. You can now rerun 
#'                                             # the model to observe any 
#'                                             # differences.
#'                                             
#' set_name(model_sirconn, 'My Epi-Model') # This function allows for setting 
#'                                         # a name for the model. Running the 
#'                                         # model object, the name of the model
#'                                         # is now reflected next to "Name of 
#'                                         # the model". 
#'                                         
#' get_name(model_sirconn) # Returns the set name of the model. 
#' 
#' get_n_viruses(model_sirconn) # Returns the number of viruses in the model. 
#'                               # In this case, there is only one virus:
#'                               # "COVID-19". 
#'                               
#' get_n_tools(model_sirconn) # Returns the number of tools in the model. In 
#'                            # this case, there are zero tools.  
#'                            
#' get_ndays(model_sirconn) # Returns the length of the simulation in days. This
#'                          # will match "ndays" within the "run" function. 
#'                          
#' get_n_replicates(model_sirconn) # Returns the number of replicates of the 
#'                                 # model. 
#'                                 
#' size(model_sirconn) # Returns the population size in the model. In this case, 
#'                     # there are 10,000 agents in the model. 
#' # Set Agents Data
#' # First, your data matrix must have the same number of rows as agents in the 
#' # model. Below is a generated matrix which will be passed into the 
#' # "set_agents_data" function. 
#' data <- matrix(data=runif(20000, min=0, max=100), nrow=10000, ncol=2)                   
#' set_agents_data(model_sirconn, data)
#' get_agents_data_ncols(model_sirconn) # Returns number of columns 
#' 
#' get_virus(model_sirconn, 0) # Returns information about the first virus in 
#'                             # the model (index begins at 0).
#' 
#' add_tool(model_sirconn, tool("Vaccine", .9, .9, .5, 1), proportion = .5)                            
#' get_tool(model_sirconn, 0) # Returns information about the first tool in the 
#'                            # model. In this case, there are no tools so an 
#'                            # error message will occur. 
queuing_on <- function(x) UseMethod("queuing_on")

#' @export
queuing_on.epiworld_sirconn <- function(x) {
  warning("SIR Connected models do not have queue.")
  invisible(x)
}

#' @export
queuing_on.epiworld_seirconn <- function(x) {
  warning("SEIR Connected models do not have queue.")
  invisible(x)
}

#' @export
queuing_on.epiworld_model <- function(x) {
  invisible(queuing_on_cpp(x))
}

#' @name epiworld-methods
#' @export
queuing_off <- function(x) UseMethod("queuing_off")

#' @export
queuing_off.epiworld_model <- function(x) {
  invisible(queuing_off_cpp(x))
}

#' @name epiworld-methods
#' @export
#' @returns 
#' - The `verbose_on` and `verbose_off` functions return the same model, however 
#' `verbose_off` returns the model with no progress bar. 
#' @details
#' The `verbose_on` and `verbose_off` functions activate and deactivate printing
#' progress on screen, respectively. Both functions return the model (`x`) invisibly.
verbose_off <- function(x) UseMethod("verbose_off")

#' @export
verbose_off.epiworld_model <- function(x) {
  invisible(verbose_off_cpp(x))
}

#' @name epiworld-methods
#' @export
verbose_on <- function(x) UseMethod("verbose_on")

#' @export
verbose_on.epiworld_model <- function(x) {
  verbose_on_cpp(x)
  invisible(x)
}

#' @export
#' @returns 
#' - The `run` function returns the simulated model of class `epiworld_model`. 
#' @rdname epiworld-methods
run <- function(model, ndays, seed = sample.int(1e4, 1)) UseMethod("run")

#' @export
run.epiworld_model <- function(model, ndays, seed = sample.int(1e4, 1)) {
  run_cpp(model, ndays, seed)
  invisible(model)
}

#' @export
print.epiworld_model <- function(x, ...) {
  print_cpp(x, lite = TRUE)
  invisible(x)
}

#' @export 
#' @returns
#' - The `summary` function prints a more detailed view of the model, and returns the same model invisibly.
#' @rdname epiworld-methods
#' @param object Object of class `epiworld_model`.
#' @param ... Additional arguments. 
summary.epiworld_model <- function(object, ...) {
  print_cpp(object, lite = FALSE)
  invisible(object)
}

#' @export
#' @returns
#' - The `get_states` function returns the unique states found in a model. 
#' @rdname epiworld-methods
get_states <- function(x) UseMethod("get_states")

#' @export
get_states.epiworld_model <- function(x) get_states_cpp(x)

#' @export
#' @param pname String. Name of the parameter.
#' @returns 
#' - The `get_param` function returns a selected parameter from the model object 
#' of class `epiworld_model`.
#' @rdname epiworld-methods
get_param <- function(x, pname) UseMethod("get_param")

#' @export
get_param.epiworld_model <- function(x, pname) {
  get_param_cpp(x, pname)
}


#' @export
#' @param pval Numeric. Value of the parameter.
#' @returns 
#' - The `set_param` function does not return a value but instead alters a
#'  parameter value. 
#' @rdname epiworld-methods
set_param <- function(x, pname, pval) UseMethod("set_param")

#' @export
set_param.epiworld_model <- function(x, pname, pval) {
  invisible(set_param_cpp(x, pname, pval))
  invisible(x)
}

#' @export
#' @param mname String. Name of the model.
#' @returns 
#' - The `set_name` function does not return a value but instead alters an object 
#' of `epiworld_model`.
#' @rdname epiworld-methods
set_name <- function(x, mname) UseMethod("set_name")

#' @export
set_name.epiworld_model <- function(x, mname) {
  set_name_cpp(x, mname)
  invisible(x)
}

#' @export
#' @returns
#' - `get_name` returns the name of the model.
#' @rdname epiworld-methods
get_name <- function(x) UseMethod("get_name")

#' @export
get_name.epiworld_model <- function(x) {
  get_name_cpp(x)
}

#' @export 
#' @rdname epiworld-methods
#' @returns
#' - `get_n_viruses` returns the number of viruses of the model.
get_n_viruses <- function(x) UseMethod("get_n_viruses")

#' @export
get_n_viruses.epiworld_model <- function(x) get_n_viruses_cpp(x)


#' @export 
#' @rdname epiworld-methods
#' @returns
#' - `get_n_tools` returns the number of tools of the model.
get_n_tools <- function(x) UseMethod("get_n_tools")

#' @export
get_n_tools.epiworld_model <- function(x) get_n_tools_cpp(x)


#' @export 
#' @rdname epiworld-methods
#' @returns
#' - `get_ndays` returns the number of days of the model.
get_ndays <- function(x) UseMethod("get_ndays")

#' @export
get_ndays.epiworld_model <- function(x) get_ndays_cpp(x)


#' @export 
#' @rdname epiworld-methods
#' @returns 
#' - `get_n_replicates` returns the number of replicates of the model.
get_n_replicates <- function(x) UseMethod("get_n_replicates")

#' @export
get_n_replicates.epiworld_model <- function(x) get_n_replicates_cpp(x)


#' @export
#' @rdname epiworld-methods
#' @returns 
#' - `size.epiworld_model` returns the number of agents in the model.
#' 
size <- function(x) UseMethod("size")

#' @export
size.epiworld_model <- function(x) size_cpp(x)


#' @export
#' @param data A numeric matrix.
#' @returns
#' - The 'set_agents_data' function returns an object of class DataFrame. 
#' @rdname epiworld-methods
set_agents_data <- function(model, data) {
  
  if (!inherits(data, "matrix") | mode(data) != "numeric")
    stop("-data- must be a numeric (mode) matrix (class).")
  
  if (size(model) != nrow(data))
    stop(
      "The number of rows in -data- (", nrow(data),
      ") doesn't match the number of agents in the model (",
      size(model), ")."
      )
  
  invisible(set_agents_data_cpp(model = model, data = data, ncols = ncol(data)))
  
}

#' @export
#' @returns 
#' - 'get_agents_data_ncols' returns the number of columns in the model dataframe. 
#' @rdname epiworld-methods
get_agents_data_ncols <- function(model) {
  
  get_agents_data_ncols_cpp(model)
  
}

#' @export
#' @param virus_pos Integer. Relative location (starting from 0) of the virus
#' in the model
#' @returns 
#' - 'get_virus' returns a [virus].
#' @rdname epiworld-methods
get_virus <- function(model, virus_pos) {
  structure(
    get_virus_model_cpp(model, virus_pos),
    class = c("epiworld_virus")
  )
}

#' @export
#' @param tool_pos Integer. Relative location (starting from 0) of the tool
#' in the model
#' @returns 
#' - `get_tool` returns a [tool].
#' @rdname epiworld-methods
get_tool <- function(model, tool_pos) {
  structure(
    get_tool_model_cpp(model, tool_pos),
    class = "epiworld_tool"
  )
}

#' @export 
#' @param proportions Numeric vector. Proportions in which agents will be
#' distributed (see details).
#' @return 
#' - `inital_states` returns the model with an updated initial state.
#' @rdname epiworld-methods
initial_states <- function(model, proportions) {

  stopifnot_model(model)
  invisible(initial_states_cpp(model, proportions))

}

#' @rdname epiworld-methods
#' @export
#' @details `epiworld_model` objects are pointers to an underlying C++ class
#' in `epiworld`. To generate a copy of a model, use `clone_model`, otherwise,
#' the assignment operator will only copy the pointer.
#' @return
#' - `clone_model` returns a copy of the model.
clone_model <- function(model) {
  stopifnot_model(model)
  structure(
    clone_model_cpp(model),
    class = class(model)
  )
}

