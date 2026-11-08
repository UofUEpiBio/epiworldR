#' Models by name: the registry
#'
#' Every built-in model whose arguments are numbers is built by
#' `make_model()`, from its id and the values of its arguments. The named
#' constructors ([ModelSIR()], [ModelSEIR()], ...) call it.
#'
#' @param id Character scalar. The model, as in `list_models()` (e.g., `"SIR"`).
#' @param values Named list (or named numeric vector) of argument values, by
#' their names in `get_model_info(id)$params$arg`. The arguments left out take
#' their defaults.
#' @param virus_name Character scalar. The name of the model's virus (`""`: the
#' default). Only for the models whose `virus_name` is not empty.
#'
#' @details
#' `make_model()` builds the model only: it has no layers, which are added with
#' [add_layer()]. A virus added to it is a strain of the model's own (see
#' [virus()]).
#'
#' `get_model_info()` lists, for each argument: its name (`arg`), the model
#' parameter it sets (`name`, empty for `n` and `prevalence`), its default, its
#' limits (checked before the model is built), whether it is a whole number,
#' and a description.
#'
#' [ModelSIRLogit()] and [ModelDiffNet()] take data and keep their own
#' constructors.
#'
#' @returns
#' - `make_model()` returns a model, of class `epiworld_<id>` and
#'   [epiworld_model].
#' - `list_models()` returns a data frame with the models' ids, labels,
#'   descriptions and default virus names.
#' - `get_model_info()` returns a list with the same fields and `params`, a data
#'   frame of the model's arguments.
#' @examples
#' list_models()
#' get_model_info("SEIR")$params
#'
#' model <- make_model(
#'   "SEIR",
#'   list(n = 1000, transmission_prob = 0.2, latent_period = 3),
#'   virus_name = "Flu"
#' )
#' add_layer(model, layer_network("contacts", graph_smallworld(8, 0.1)))
#' run(model, ndays = 50, seed = 1)
#' model
#' @export
#' @concept general-models
make_model <- function(id, values = list(), virus_name = "") {

  stopifnot_string(id)
  stopifnot_string(virus_name)

  if (length(values) && (is.null(names(values)) || any(names(values) == "")))
    stop("-values- must be named, by argument (see get_model_info()).", call. = FALSE)

  for (arg in names(values)) {
    v <- values[[arg]]
    if (!is.numeric(v) || length(v) != 1L || is.na(v))
      stop(sprintf("-%s- must be a number (not NA).", arg), call. = FALSE)
  }
  values <- unlist(values)

  # as.double() drops the names
  v <- as.double(values)
  names(v) <- names(values)

  structure(
    make_model_cpp(id, v, virus_name),
    class = c(paste0("epiworld_", tolower(id)), "epiworld_model")
  )

}

#' @export
#' @rdname make_model
list_models <- function() {
  list_models_cpp()
}

#' @export
#' @rdname make_model
get_model_info <- function(id) {
  stopifnot_string(id)
  info <- get_model_info_cpp(id)
  info$params <- as.data.frame(info$params)
  info
}
