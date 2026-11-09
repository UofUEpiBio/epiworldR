#' Defunct functions in epiworldR
#'
#' These functions were removed in epiworldR 1.0.0, when the population and
#' the contacts of a model moved to its constructor and its contact layers
#' (see [contact-layers]), and the course of infection moved to the virus's
#' natural history (its stages and transitions, see [virus()]). A model's
#' layers describe who meets whom, so the models that differed only in where
#' contacts came from are gone too; and one vocabulary ([assign]) places
#' viruses, tools, initial stages and groups.
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
  "(one group per agent) or a placement (e.g., assign_counts(sizes));",
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

# The natural history (epiworld NH-1, NH-2) replaced the states and their
# update functions; viruses carry their stages and transitions.
defunct_state <- paste(
  "States are gone: the virus carries the stages of the disease and the",
  "transitions between them (add_stage(), add_infection(), add_rate(),",
  "add_after(); see help(\"virus\"))."
)

# Placement (decision 14 of natural-history.md): one vocabulary, assign_*()
defunct_assign_virus <- paste(
  "Place the virus when adding it: add_virus(model, virus, assign), with",
  "assign_share(p), assign_count(n), assign_agents(ids) or",
  "assign_group_shares(layer, p) (see help(\"assign\"))."
)

defunct_assign_tool <- paste(
  "Place the tool when adding it: add_tool(model, tool, assign), with",
  "assign_share(p), assign_count(n), assign_agents(ids) or",
  "assign_group_shares(layer, p) (see help(\"assign\"))."
)

defunct_virus_value <- function(what) {
  paste(
    "A virus holds no values: its transitions name model parameters.",
    sprintf("Declare the %s as a transition", what),
    "(add_infection(), add_rate() or add_after(), with a number, a parameter",
    "name, value_fun() or rate_fun_logit()); a strain reads other parameters",
    "with use_param()."
  )
}

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
  distribute_tool_to_entities = defunct_assign_tool,
  distribute_virus_to_entities = defunct_assign_virus,
  add_state = defunct_state,
  set_state_function = defunct_state,
  update_fun_rate = defunct_state,
  update_fun_susceptible = defunct_state,
  virus_set_state = defunct_state,
  set_infectious_states = paste(
    "A stage transmits if it is declared infectious:",
    "add_stage(virus, label, infectious = TRUE)."
  ),
  get_infectious_states = "Use get_stages(virus) (its -infectious- column).",
  get_states = "Use get_stages(model) and get_statuses(model).",
  get_state = "Use get_stage(model, agent) and get_status(model, agent).",
  get_agents_states = "Use get_agents_stages(model) and get_agents_statuses(model).",
  change_state = paste(
    "An agent moves only along a declared transition: fire(model, transition,",
    "agent_id) (e.g., an external one, add_external())."
  ),
  initial_states = paste(
    "Use set_initial(model, assign_shares(c(<stage> = <share>, ...))): the",
    "agents start in stages without the virus."
  ),
  set_prob_infecting = defunct_virus_value("probability of infection"),
  set_prob_infecting_ptr = defunct_virus_value("probability of infection"),
  set_prob_infecting_fun = defunct_virus_value("probability of infection"),
  set_prob_recovery = defunct_virus_value("recovery"),
  set_prob_recovery_ptr = defunct_virus_value("recovery"),
  set_prob_recovery_fun = defunct_virus_value("recovery"),
  set_prob_death = defunct_virus_value("death"),
  set_prob_death_ptr = defunct_virus_value("death"),
  set_prob_death_fun = defunct_virus_value("death"),
  set_incubation = defunct_virus_value("incubation"),
  set_incubation_ptr = defunct_virus_value("incubation"),
  set_incubation_fun = defunct_virus_value("incubation"),
  virus_fun_logit = "Use rate_fun_logit(vars, coefs, model) as a transition's value.",
  set_distribution_virus = defunct_assign_virus,
  distribute_virus_randomly = defunct_assign_virus,
  distribute_virus_to_set = defunct_assign_virus,
  distribute_virus_set = defunct_assign_virus,
  distribute_virus_to_groups = defunct_assign_virus,
  add_virus_n = defunct_assign_virus,
  set_distribution_tool = defunct_assign_tool,
  distribute_tool_randomly = defunct_assign_tool,
  distribute_tool_to_set = defunct_assign_tool,
  distribute_tool_to_groups = defunct_assign_tool,
  add_tool_n = defunct_assign_tool,
  grouping_random = paste(
    "Give layer_mixing() a placement whose labels are the groups:",
    "assign_counts(sizes) or assign_shares(proportions) (see help(\"assign\"))."
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

#' @export
#' @rdname epiworldR-defunct
add_state <- function(...) defunct("add_state")

#' @export
#' @rdname epiworldR-defunct
set_state_function <- function(...) defunct("set_state_function")

#' @export
#' @rdname epiworldR-defunct
update_fun_rate <- function(...) defunct("update_fun_rate")

#' @export
#' @rdname epiworldR-defunct
update_fun_susceptible <- function(...) defunct("update_fun_susceptible")

#' @export
#' @rdname epiworldR-defunct
virus_set_state <- function(...) defunct("virus_set_state")

#' @export
#' @rdname epiworldR-defunct
set_infectious_states <- function(...) defunct("set_infectious_states")

#' @export
#' @rdname epiworldR-defunct
get_infectious_states <- function(...) defunct("get_infectious_states")

#' @export
#' @rdname epiworldR-defunct
get_states <- function(...) defunct("get_states")

#' @export
#' @rdname epiworldR-defunct
get_state <- function(...) defunct("get_state")

#' @export
#' @rdname epiworldR-defunct
get_agents_states <- function(...) defunct("get_agents_states")

#' @export
#' @rdname epiworldR-defunct
change_state <- function(...) defunct("change_state")

#' @export
#' @rdname epiworldR-defunct
initial_states <- function(...) defunct("initial_states")

#' @export
#' @rdname epiworldR-defunct
set_prob_infecting <- function(...) defunct("set_prob_infecting")

#' @export
#' @rdname epiworldR-defunct
set_prob_infecting_ptr <- function(...) defunct("set_prob_infecting_ptr")

#' @export
#' @rdname epiworldR-defunct
set_prob_infecting_fun <- function(...) defunct("set_prob_infecting_fun")

#' @export
#' @rdname epiworldR-defunct
set_prob_recovery <- function(...) defunct("set_prob_recovery")

#' @export
#' @rdname epiworldR-defunct
set_prob_recovery_ptr <- function(...) defunct("set_prob_recovery_ptr")

#' @export
#' @rdname epiworldR-defunct
set_prob_recovery_fun <- function(...) defunct("set_prob_recovery_fun")

#' @export
#' @rdname epiworldR-defunct
set_prob_death <- function(...) defunct("set_prob_death")

#' @export
#' @rdname epiworldR-defunct
set_prob_death_ptr <- function(...) defunct("set_prob_death_ptr")

#' @export
#' @rdname epiworldR-defunct
set_prob_death_fun <- function(...) defunct("set_prob_death_fun")

#' @export
#' @rdname epiworldR-defunct
set_incubation <- function(...) defunct("set_incubation")

#' @export
#' @rdname epiworldR-defunct
set_incubation_ptr <- function(...) defunct("set_incubation_ptr")

#' @export
#' @rdname epiworldR-defunct
set_incubation_fun <- function(...) defunct("set_incubation_fun")

#' @export
#' @rdname epiworldR-defunct
virus_fun_logit <- function(...) defunct("virus_fun_logit")

#' @export
#' @rdname epiworldR-defunct
set_distribution_virus <- function(...) defunct("set_distribution_virus")

#' @export
#' @rdname epiworldR-defunct
distribute_virus_randomly <- function(...) defunct("distribute_virus_randomly")

#' @export
#' @rdname epiworldR-defunct
distribute_virus_to_set <- function(...) defunct("distribute_virus_to_set")

#' @export
#' @rdname epiworldR-defunct
distribute_virus_set <- function(...) defunct("distribute_virus_set")

#' @export
#' @rdname epiworldR-defunct
distribute_virus_to_groups <- function(...) defunct("distribute_virus_to_groups")

#' @export
#' @rdname epiworldR-defunct
add_virus_n <- function(...) defunct("add_virus_n")

#' @export
#' @rdname epiworldR-defunct
set_distribution_tool <- function(...) defunct("set_distribution_tool")

#' @export
#' @rdname epiworldR-defunct
distribute_tool_randomly <- function(...) defunct("distribute_tool_randomly")

#' @export
#' @rdname epiworldR-defunct
distribute_tool_to_set <- function(...) defunct("distribute_tool_to_set")

#' @export
#' @rdname epiworldR-defunct
distribute_tool_to_groups <- function(...) defunct("distribute_tool_to_groups")

#' @export
#' @rdname epiworldR-defunct
add_tool_n <- function(...) defunct("add_tool_n")

#' @export
#' @rdname epiworldR-defunct
grouping_random <- function(...) defunct("grouping_random")
