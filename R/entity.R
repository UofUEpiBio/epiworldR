
stopifnot_entity <- function(entity) {
  if (!inherits(entity, "epiworld_entity")) {
    stop("Argument 'entity' must be an entity object.")
  }
}

#' Get entities
#' 
#' Entities in `epiworld` are objects that can contain agents.
#' @param model Model object of class `epiworld_model`.
#' 
#' @details 
#' Epiworld entities are especially useful for mixing models, particularly
#' [ModelSIRMixing] and [ModelSEIRMixing].
#' 
#' @name entities
#' @export 
#' @examples 
#' # Creating a mixing model
#' mymodel <- ModelSIRMixing(
#'  name = "My model",
#'  n = 10000,
#'  prevalence = .001,
#'  contact_rate = 10,
#'  transmission_rate = .1,
#'  recovery_rate = 1/7,
#'  contact_matrix = matrix(c(.9, .1, .1, .9), 2, 2)
#' )
#' 
#' ent1 <- entity("First")
#' ent2 <- entity("Second")
#' 
#' add_entity_n(mymodel, ent1, 5000)
#' add_entity_n(mymodel, ent2, 5000)
#' 
#' run(mymodel, ndays = 100, seed = 1912)
get_entities <- function(model) {
  stopifnot_model(model)
  structure(
    get_entities_cpp(model),
    class = c("epiworld_entities"),
    model = model
  )
}

#' @export 
#' @rdname entities
#' @param x Object of class `epiworld_entities`.
#' @param i Integer index.
`[.epiworld_entities` <- function(x, i) {

  stopifnot_entity(x)
  
  if (i > get_entity_size(x)) {
    stop("Index out of bounds.")
  }

  structure(
    get_entity_cpp(x, i),
    class = c("epiworld_entity"),
    model = x$model
  )

}

#' @export 
#' @param name Character scalar. Name of the entity.
#' @return 
#' - The function `entity` creates an entity object.
#' @rdname entities
entity <- function(name) {

  structure(
    entity_cpp(name),
    class = "epiworld_entity"
  )

}

#' @export 
#' @rdname entities
#' @param entity Entity object of class `epiworld_entity`.
#' @return
#' - The function `get_entity_size` returns the number of agents in the entity.
get_entity_size <- function(entity) {
  stopifnot_entity(entity)
  get_entity_size_cpp(entity)
}

#' @export
#' @rdname entities
#' @return 
#' - The function `get_entity_name` returns the name of the entity.
get_entity_name <- function(entity) {
  stopifnot_entity(entity)
  get_entity_name_cpp(entity)
}

#' @export
#' @rdname entities
#' @param agent Agent object of class `epiworld_agent`.
#' @return 
#' - The function `entity_add_agent` adds an agent to the entity.
entity_add_agent <- function(
  entity,
  agent,
  model = attr(entity, "model")
  ) {

  stopifnot_entity(entity)
  stopifnot_agent(agent)
  entity_add_agent_cpp(entity, agent, model)

  invisible(entity)

}

#' @export
#' @rdname entities
#' @return 
#' - The function `entity_rm_agent` removes an agent from the entity.
entity_rm_agent <- function(
  entity,
  id
) {

  stopifnot_entity(entity)
  stopifnot_agent(agent)
  entity_rm_agent_cpp(entity, agent)

  invisible(entity)

}

# #' @export
# #' @rdname entities
# #' @return 
# #' - The function `rm_entity` removes an entity from the model.
# rm_entity <- function(model, id) {
  
#   stopifnot_model(model)
#   rm_entity_cpp(model, entity)

#   invisible(model)
# }

#' @export
#' @rdname entities
#' @param n Integer scalar. Number of randomly assigned agents.
add_entity_n <- function(
  model,
  entity,
  n
) {

  stopifnot_model(model)
  stopifnot_entity(entity)
  add_entity_n_cpp(model, entity, n)

  invisible(model)

}
