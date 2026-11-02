#' Defunct functions in epiworldR
#'
#' These functions were removed in epiworldR 1.0.0, when the population and
#' the contacts of a model moved to its constructor and its contact layers
#' (see [contact-layers]). A model's states describe the course of infection and
#' its layers describe who meets whom, so the models that differed only in
#' where contacts came from are gone too.
#' Calling them is an error that says what to use instead:
#'
#' `r defunct_doc()`
#'
#' @param ... Ignored.
#' @name epiworldR-defunct
#' @keywords internal
NULL

# Entities (removed in epiworld's phase 3d-3a) became the groups of a mixing
# layer.
defunct_entity <- paste(
  "Entities are gone: groups live in a mixing layer.",
  "add_layer(model, layer_mixing(name, groups, ...)), with groups a vector",
  "(one group per agent) or grouping_random(sizes, as_proportion);",
  "members_of(layer, group) and get_group(layer, agent) read them."
)

defunct_entity_fns <- c(
  "get_entities",
  "entity",
  "get_entity_size",
  "get_entity_name",
  "entity_add_agent",
  "add_entity",
  "rm_entity",
  "load_agents_entities_ties",
  "entity_get_agents",
  "distribute_entity_randomly",
  "distribute_entity_to_set",
  "set_distribution_entity",
  "add_entities_from_dataframe"
)

# What replaces each defunct function. This one table feeds the stubs and
# their help page.
defunct_table <- c(
  agents_smallworld = paste(
    "Give the population size to the model constructor (e.g.,",
    "ModelSIR(name, n, ...)) and add a network layer:",
    "add_layer(model, layer_network(\"network\", graph_smallworld(k, p)))."
  ),
  agents_from_edgelist = paste(
    "Give the population size to the model constructor and add a fixed",
    "network layer: add_layer(model, layer_network(\"network\",",
    "graph_edgelist(source, target, size, directed)))."
  ),
  agents_sbm = paste(
    "Give the population size to the model constructor and add a network",
    "layer: add_layer(model, layer_network(\"network\",",
    "graph_sbm(block_sizes, mixing_matrix)))."
  ),
  get_network = "Use get_edgelist(get_layer(model, name)).",
  ModelSIRCONN = paste(
    "Use ModelSIR() with a random-mixing layer: everyone in one group, a",
    "contact matrix of 1, and set_contact_rate(layer, \"Contact rate\")",
    "(see the example in help(\"mixing-layers\"))."
  ),
  ModelSEIRCONN = paste(
    "Use ModelSEIR() with a random-mixing layer: everyone in one group, a",
    "contact matrix of 1, and set_contact_rate(layer, \"Contact rate\")",
    "(see the example in help(\"mixing-layers\"))."
  ),
  ModelSIRDCONN = paste(
    "Use ModelSIRD() with a random-mixing layer: everyone in one group, a",
    "contact matrix of 1, and set_contact_rate(layer, \"Contact rate\")",
    "(see the example in help(\"mixing-layers\"))."
  ),
  ModelSEIRDCONN = paste(
    "Use ModelSEIRD() with a random-mixing layer: everyone in one group, a",
    "contact matrix of 1, and set_contact_rate(layer, \"Contact rate\")",
    "(see the example in help(\"mixing-layers\"))."
  ),
  ModelSIRMixing = paste(
    "Use ModelSIR() and add_layer(model, layer_mixing(name, groups,",
    "contact_matrix = contact_matrix))."
  ),
  ModelSEIRMixing = paste(
    "Use ModelSEIR() and add_layer(model, layer_mixing(name, groups,",
    "contact_matrix = contact_matrix))."
  ),
  ModelSEIRMixingQuarantine = paste(
    "Use ModelSEIRQuarantine() and add_layer(model, layer_mixing(name, groups,",
    "contact_matrix = contact_matrix)); it traces every layer it is given."
  ),
  setNames(rep(defunct_entity, length(defunct_entity_fns)), defunct_entity_fns),
  distribute_tool_to_entities = paste(
    "Use distribute_tool_to_groups(layer, prevalence, as_proportion) with the",
    "name of a mixing layer."
  ),
  distribute_virus_to_entities = paste(
    "Use distribute_virus_to_groups(layer, prevalence, as_proportion) with the",
    "name of a mixing layer."
  )
)

defunct <- function(old) {
  .Defunct(msg = sprintf(
    "'%s()' is defunct as of epiworldR 1.0.0. %s See help(\"epiworldR-defunct\").",
    old, defunct_table[[old]]
  ))
}

defunct_doc <- function() {
  paste0(
    "- `", names(defunct_table), "()`: ", defunct_table,
    collapse = "\n"
  )
}

#' @export
#' @rdname epiworldR-defunct
agents_smallworld <- function(...) defunct("agents_smallworld")

#' @export
#' @rdname epiworldR-defunct
agents_from_edgelist <- function(...) defunct("agents_from_edgelist")

#' @export
#' @rdname epiworldR-defunct
agents_sbm <- function(...) defunct("agents_sbm")

#' @export
#' @rdname epiworldR-defunct
get_network <- function(...) defunct("get_network")

#' @export
#' @rdname epiworldR-defunct
ModelSIRCONN <- function(...) defunct("ModelSIRCONN")

#' @export
#' @rdname epiworldR-defunct
ModelSEIRCONN <- function(...) defunct("ModelSEIRCONN")

#' @export
#' @rdname epiworldR-defunct
ModelSIRDCONN <- function(...) defunct("ModelSIRDCONN")

#' @export
#' @rdname epiworldR-defunct
ModelSEIRDCONN <- function(...) defunct("ModelSEIRDCONN")

#' @export
#' @rdname epiworldR-defunct
ModelSIRMixing <- function(...) defunct("ModelSIRMixing")

#' @export
#' @rdname epiworldR-defunct
ModelSEIRMixing <- function(...) defunct("ModelSEIRMixing")

#' @export
#' @rdname epiworldR-defunct
ModelSEIRMixingQuarantine <- function(...) defunct("ModelSEIRMixingQuarantine")

#' @export
#' @rdname epiworldR-defunct
get_entities <- function(...) defunct("get_entities")

#' @export
#' @rdname epiworldR-defunct
entity <- function(...) defunct("entity")

#' @export
#' @rdname epiworldR-defunct
get_entity_size <- function(...) defunct("get_entity_size")

#' @export
#' @rdname epiworldR-defunct
get_entity_name <- function(...) defunct("get_entity_name")

#' @export
#' @rdname epiworldR-defunct
entity_add_agent <- function(...) defunct("entity_add_agent")

#' @export
#' @rdname epiworldR-defunct
add_entity <- function(...) defunct("add_entity")

#' @export
#' @rdname epiworldR-defunct
rm_entity <- function(...) defunct("rm_entity")

#' @export
#' @rdname epiworldR-defunct
load_agents_entities_ties <- function(...) defunct("load_agents_entities_ties")

#' @export
#' @rdname epiworldR-defunct
entity_get_agents <- function(...) defunct("entity_get_agents")

#' @export
#' @rdname epiworldR-defunct
distribute_entity_randomly <- function(...) defunct("distribute_entity_randomly")

#' @export
#' @rdname epiworldR-defunct
distribute_entity_to_set <- function(...) defunct("distribute_entity_to_set")

#' @export
#' @rdname epiworldR-defunct
set_distribution_entity <- function(...) defunct("set_distribution_entity")

#' @export
#' @rdname epiworldR-defunct
add_entities_from_dataframe <- function(...) defunct("add_entities_from_dataframe")

#' @export
#' @rdname epiworldR-defunct
distribute_tool_to_entities <- function(...) defunct("distribute_tool_to_entities")

#' @export
#' @rdname epiworldR-defunct
distribute_virus_to_entities <- function(...) defunct("distribute_virus_to_entities")
