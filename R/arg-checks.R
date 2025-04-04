# This file contains functions for checking the arguments
# of the exported package functions.


# Checks if argument contains NA values
stopifany_na <- function(x) {
  if (any(is.na(x))) {
    stop(paste(match.call()$x, "must not contain NA values."))
  }
}

# Checks if argument is a string
stopifnot_string <- function(x) {
  if (!is.character(x)) {
    stop(
      paste(match.call()$x, "must be a string, but is of class(es): "),
      paste(class(x), collapse = ", ")
    )
  }
}

# Checks if argument is an integer
stopifnot_int <- function(x) {
  if (!is.numeric(x) || !all.equal(x, as.integer(x))) {
    stop(
      paste(match.call()$x, "must be an integer, but is of class(es): "),
      paste(class(x), collapse = ", ")
    )
  }
}

# Checks if argument is a double
stopifnot_double <- function(x) {
  if (!is.numeric(x)) {
    stop(
      paste(match.call()$x, "must be a double, but is of class(es): "),
      paste(class(x), collapse = ", ")
    )
  }
}

# Checks if argument is a boolean
stopifnot_bool <- function(x) {
  if (is.na(x)) {
    stop(paste(match.call()$x, "must not be NA."))
  }

  if (!is.logical(x)) {
    stop(
      paste(match.call()$x, "must be a boolean, but is of class(es): "),
      paste(class(x), collapse = ", ")
    )
  }
}

# Checks if argument is a vector of numeric values
stopifnot_numvector <- function(x) {
  if (!is.numeric(x) || !is.vector(x)) {
    stop(
      paste(match.call()$x, "must be a numeric vector, but is of class(es): "),
      paste(class(x), collapse = ", ")
    )
  }
}

# Checks if argument is a vector of string values
stopifnot_stringvector <- function(x) {
  if (!is.character(x) || !is.vector(x)) {
    stop(
      paste(match.call()$x, "must be a string vector, but is of class(es): "),
      paste(class(x), collapse = ", ")
    )
  }
}


# Checks if model object is of class "epiworld_model"
stopifnot_model <- function(model) {
  if (!inherits(model, "epiworld_model")) {
    stop(
      "The -model- object must be of class 'epiworld_model'. ",
      "The object passed to the function is of class(es): ",
      paste(class(model), collapse = ", ")
    )
  }
}

# Checks if agent object is of class "epiworld_agent"
stopifnot_agent <- function(agent) {
  if (!inherits(agent, "epiworld_agent"))
    stop(
      "The -agent- object must be of class 'epiworld_agent'. ",
      "The object passed to the function is of class(es): ",
      paste(class(agent), collapse = ", ")
    )
}

# checks if entity object is of class "epiworld_entity"
stopifnot_entity <- function(entity) {
  if (!inherits(entity, "epiworld_entity")) {
    stop(
      "The -entity- object must be of class 'epiworld_entity'. ",
      "The object passed to the function is of class(es): ",
      paste(class(entity), collapse = ", ")
    )
  }
}

# Checks if distfun object is of class "epiworld_distribution_entity"
stopifnot_entity_distfun <- function(distfun) {
  if (!inherits(distfun, "epiworld_distribution_entity")) {
    stop("Argument 'distfun' must be a distribution function.")
  }
}

# Checks if object is of class "epiworld_lfmcmc"
stopifnot_lfmcmc <- function(x) {
  # Catching the value of x
  nam <- match.call()$x

  if (!inherits(x, "epiworld_lfmcmc"))
    stop(nam, " must be an object of class epiworld_lfmcmc")

}

# Checks if tool object is of class "epiworld_tool"
stopifnot_tool <- function(tool) {
  if (!inherits(tool, "epiworld_tool")) {
    stop(
      "The -tool- object must be of class 'epiworld_tool'. ",
      "The object passed to the function is of class(es): ",
      paste(class(tool), collapse = ", ")
    )
  }
}

# Checks if tfun is of class "epiworld_tool_fun"
stopifnot_tfun <- function(tfun) {
  if (!inherits(tfun, "epiworld_tool_fun")) {
    stop(
      "The -tfun- object must be of class 'epiworld_tool_fun'. ",
      "The object passed to the function is of class(es): ",
      paste(class(tfun), collapse = ", ")
    )
  }
}

# Checks if tool_distfun is of class "epiworld_tool_distfun"
stopifnot_tool_distfun <- function(tool_distfun) {
  if (!inherits(tool_distfun, "epiworld_tool_distfun")) {
    stop(
      "The -tool_distfun- object must be of class 'epiworld_tool_distfun'. ",
      "The object passed to the function is of class(es): ",
      paste(class(tool_distfun), collapse = ", ")
    )
  }
}

# Checks if virus object is of class "epiworld_virus"
stopifnot_virus <- function(virus) {
  if (!inherits(virus, "epiworld_virus")) {
    stop(
      "The -virus- object must be of class 'epiworld_virus'. ",
      "The object passed to the function is of class(es): ",
      paste(class(virus), collapse = ", ")
    )
  }
}

# Checks if vfun is of class "epiworld_virus_fun"
stopifnot_vfun <- function(vfun) {
  if (!inherits(vfun, "epiworld_virus_fun")) {
    stop(
      "The -vfun- object must be of class 'epiworld_virus_fun'. ",
      "The object passed to the function is of class(es): ",
      paste(class(vfun), collapse = ", ")
    )
  }
}

# Checks if virus_distfun is of class "epiworld_virus_distfun"
stopifnot_virus_distfun <- function(virus_distfun) {
  if (!inherits(virus_distfun, "epiworld_virus_distfun")) {
    stop(
      "The -virus_distfun- object must be of class 'epiworld_virus_distfun'. ",
      "The object passed to the function is of class(es): ",
      paste(class(virus_distfun), collapse = ", ")
    )
  }
}

# Checks if argument is of class "epiworld_modeldiagram"
stopifnot_modeldiagram <- function(x) {
  if (!inherits(x, "epiworld_modeldiagram")) {
    stop(
      paste(match.call()$x, "must be an object of class 'epiworld_modeldiagram', but is of class(es): "),
      paste(class(x), collapse = ", ")
    )
  }
}
