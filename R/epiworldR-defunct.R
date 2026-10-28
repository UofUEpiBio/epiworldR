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
