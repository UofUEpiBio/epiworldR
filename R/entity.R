
stopifnot_entity <- function(entity) {
  if (!inherits(entity, "epiworld_entity")) {
    stop("Argument 'entity' must be an entity object.")
  }
}

stopifnot_entity_distfun <- function(distfun) {
  if (!inherits(distfun, "epiworld_distribution_entity")) {
    stop("Argument 'distfun' must be a distribution function.")
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
#' ent1 <- entity("First", 5000, FALSE)
#' ent2 <- entity("Second", 5000, FALSE)
#' 
#' mymodel |>
#'   add_entity(ent1) |>
#'   add_entity(ent2)
#' 
#' run(mymodel, ndays = 100, seed = 1912)
#' 
#' summary(mymodel)
get_entities <- function(model) {

  stopifnot_model(model)
  structure(
    lapply(
      get_entities_cpp(model), \(e) {
        structure(
          e,
          class = c("epiworld_entity"),
          model = model
        )
      }
    ),
    class = c("epiworld_entities")
  )
}

#' @export
print.epiworld_entities <- function(x, ...) {
  cat("A collection of ", length(x), " entities.\n")
  invisible(x)
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
#' @param prevalence Numeric scalar. Prevalence of the entity.
#' @param as_proportion Logical scalar. If `TRUE`, `prevalence` is interpreted
#' as a proportion.
#' @param to_unassigned Logical scalar. If `TRUE`, the entity is added to the
#' unassigned pool.
#' @return 
#' - The function `entity` creates an entity object.
#' @rdname entities
entity <- function(name, prevalence, as_proportion, to_unassigned = TRUE) {

  structure(
    entity_cpp(
      name,
      as.double(prevalence),
      as.logical(as_proportion),
      as.logical(to_unassigned)
      ),
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
#' @param agent_id Integer scalar. Agent id to remove from the entity.
#' @return 
#' - The function `entity_rm_agent` removes an agent from the entity.
entity_rm_agent <- function(
  entity,
  agent_id
) {

  stopifnot_entity(entity)
  entity_rm_agent_cpp(entity, agent_id)

  invisible(entity)

}

#' @export
#' @rdname entities
#' @param id Integer scalar. Entity id to remove (starting from zero).
#' @return 
#' - The function `rm_entity` removes an entity from the model.
rm_entity <- function(model, id) {
  
  stopifnot_model(model)
  rm_entity_cpp(model, entity)

  invisible(model)
}

#' @export
#' @rdname entities
add_entity <- function(
  model,
  entity
) {

  stopifnot_model(model)
  stopifnot_entity(entity)
  add_entity_cpp(
    model,
    entity
    )

  invisible(model)

}

#' @export 
#' @rdname entities
#' @param agents_id Integer vector. 
#' @param entities_id Integer vector. 
#' @return 
#' - The function `load_agents_entities_ties` loads agents into entities.
load_agents_entities_ties <- function(
  model,
  agents_id,
  entities_id
) {

  stopifnot_model(model)
  if (!inherits(agents_id, "integer")) {
    stop("Argument 'agents_id' must be an integer.")
  }

  if (!inherits(entities_id, "integer")) {
    stop("Argument 'entities_id' must be an integer.")
  }

  load_agents_entities_ties_cpp(model, agents_id, entities_id)

  invisible(model)

}

#' @export 
#' @rdname entities
#' @return 
#' - The function `entity_get_agents` returns an integer vector with the agents
#'   in the entity (ids).
entity_get_agents <- function(entity) {

  stopifnot_entity(entity)
  entity_get_agents_cpp(entity)

}

#' @export 
print.epiworld_entity <- function(x, ...) {
  print_entity_cpp(x)
  invisible(x)
}

#' @export
#' @param prevalence Numeric scalar. Prevalence of the entity.
#' @param as_proportion Logical scalar. If `TRUE`, `prevalence` is interpreted
#' as a proportion.
#' @rdname entities
distribute_entity_randomly <- function(
  prevalence,
  as_proportion,
  to_unassigned = TRUE
) {

  structure(
      distribute_entity_randomly_cpp(
      as.double(prevalence),
      as.logical(as_proportion),
      as.logical(to_unassigned)
    ),
    class = "epiworld_distribution_entity"
  )

}

#' @export
#' @param agents_ids Integer vector. Ids of the agents to distribute.
#' @rdname entities
distribute_entity_to_set <- function(
  agents_ids
) {

  structure(
    distribute_entity_a_set_cpp(
      as.integer(agents_ids)
    ),
    class = "epiworld_distribution_entity"
  )

}

#' @export 
#' @rdname entities
#' @param distfun Distribution function object of class `epiworld_distribution_entity`.
set_distribution_entity <- function(
  entity,
  distfun
) {

  stopifnot_entity(entity)
  stopifnot_entity_distfun(distfun)
  set_distribution_entity_cpp(entity, distribution)

  invisible(entity)

}