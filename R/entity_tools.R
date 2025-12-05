#' Add entities to a model according to a data.frame
#'
#' Helper function that facilities creating entities and adding them to
#' models. It is a wrapper of [add_entity()].
#'
#' @param model An [epiworld_model] object.
#' @param entities A `data.frame` with the entities to add. It must contain two
#' columns: entity names (character) and size (proportion or integer).
#' @param col_name The name of the column in `entities` that contains the
#' entity names.
#' @param col_number The name of the column in `entities` that contains the
#' entity sizes (either as proportions or integers).
#' @param ... Further arguments passed to [add_entity()]
#' @returns
#' Inivisible the model with the added entities.
#' @examples
#' # Start off creating three entities.
#' # Individuals will be distributed randomly between the three.
#' entities <- data.frame(
#'   name = c("Pop 1", "Pop 2", "Pop 3"),
#'   size = rep(3e3, 3)
#' )
#'
#' # Row-stochastic matrix (rowsums 1)
#' cmatrix <- c(
#'   c(0.9, 0.05, 0.05),
#'   c(0.1, 0.8, 0.1),
#'   c(0.1, 0.2, 0.7)
#' ) |> matrix(byrow = TRUE, nrow = 3)
#'
#' flu_model <- ModelSIRMixing(
#'   name              = "Flu",
#'   n                 = 9e3,
#'   prevalence        = 1 / 9e3,
#'   contact_rate      = 20,
#'   transmission_rate = 0.1,
#'   recovery_rate     = 1 / 7,
#'   contact_matrix    = cmatrix
#' ) |>
#'   add_entities_from_dataframe(
#'     entities = entities,
#'     col_name = "name",
#'     col_number = "size",
#'     # This is passed to `add_entity()`
#'     as_proportion = FALSE
#'   )
#' @export
add_entities_from_dataframe <- function(
  model,
  entities,
  col_name,
  col_number,
  ...
) {
  # Basic checker
  stopifnot_model(model)
  stopifnot_columns(entities, c(col_name, col_number))
  stopifnot_int(entities[[col_number]], lb = 1L)

  # Creating sequences of agent ids from 0 to n-1
  # So we can distribute them to entities
  roll_sum <- 0L

  # Iterating through the rows
  for (i in seq_len(nrow(entities))) {
    # Setting the distribution set
    from <- roll_sum
    to   <- roll_sum + entities[[col_number]][i] - 1L

    e <- entity(
      name = entities[[col_name]][i],
      prevalence = entities[[col_number]][i],
      ...
    )

    # Setting the distribution to the set
    set_distribution_entity(
      entity = e,
      distfun = distribute_entity_to_set(from:to)
    )

    add_entity(model, e)

    roll_sum <- to + 1L

  }

  invisible(model)

}
