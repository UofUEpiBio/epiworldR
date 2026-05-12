#' Visualize results from `run_multiple()`
#' @param x An object of class `epiworld_multiple_save`
#' @param y Ignored.
#' @param ... Further arguments passed to the method.
#' @details
#' In the case of `transition` objects, the function will collapse the
#' transition matrix and call [draw_mermaid_from_matrix()] followed by
#' `plot()`.
#' @export
#' @concept plot-epidemic-curves
plot.epiworld_multiple_save <- function(x, y = NULL, ...) {
  # what <- attr(x, "what")
  invisible(lapply(x, plot, ...))

}

#' Compute time-series confidence intervals area
#'
#' This function computes the area for confidence intervals of time-series
#' data, which can be used for visualizing uncertainty in time-series plots.
#' It takes a numeric vector `x`, a grouping vector `group` that indicates
#' which time point each value in `x` corresponds to, and an `alpha` level for
#' the confidence interval. The function returns a list containing two data
#' frames: one for the area of the confidence interval (suitable for polygon
#' plotting) and one for the median line of the time-series.
#'
#' @param x Numeric vector.
#' @param group Grouping vector of the same length as `x`.
#' @param alpha Numeric scalar. Size of the two-sided interval tail.
#' @return A list with two data frames: `area` (for polygons) and `line`
#' (for the median time-series).
#' @examples
#' # Simulating random walks
#' set.seed(123)
#' dat <- lapply(1:200, \(i) {
#'   data.frame(
#'     x = cumsum(rnorm(10)),
#'     time = 1:10
#'   )
#' }) |> do.call(what = "rbind")
#' ans <- compute_ts_ci_area(dat$x, dat$time)
#' plot(ans$line, type = "b", ylim = c(-6, 6))
#' polygon(x = ans$area$x, y = ans$area$y)
#' @export
#' @concept misc-functions
compute_ts_ci_area <- function(
  x,
  group,
  alpha = .05
) {

  alpha_half <- alpha / 2

  res <- tapply(x, list(group), FUN = \(x) {
    stats::quantile(
      x, prob = c(
        alpha_half,
        .5,
        1 - alpha_half
      )
    )
  })

  res <- cbind(
    data.frame(index = as.integer(names(res))),
    do.call(rbind, res)
  )

  # Ensuring data is sorted + naming
  res <- res[order(res[, 1]), ]

  # Now, creating curve points
  list(
    area = data.frame(
      x = c(res[[1]], rev(res[[1]])),
      y = c(res[[4]], rev(res[[2]]))
    ),
    line = data.frame(
      x = res[[1]],
      y = res[[3]]
    )
  )
}

#' Plot transition matrices from `run_multiple_get_results()`
#'
#' This function takes the transition data frame from
#' `run_multiple_get_results()`, particularly, the `transition` component,
#' and visualizes it as a transition matrix using `draw_mermaid_from_matrix()`.
#' @param x A data frame in transition format, usually
#' `run_multiple_get_results(model)$transition`.
#' @param ... Further arguments passed to [draw_mermaid_from_matrix()].
#' @details
#' The function is automatically called when plotting an object of class
#' `epiworld_multiple_save`. It constructs a transition matrix from the
#' provided data frame, normalizes the rows to represent probabilities, and
#' then visualizes it using a mermaid diagram. The resulting plot illustrates
#' the transitions between different states in the model, with the thickness of
#' the arrows corresponding to the transition.
#' @export
#' @concept epiworld-model-diagrams
plot_multiple_transition <- function(x, ...) {
  # Figuring out the states
  states <- unique(c(x$from, x$to)) |>
    sort()

  # Creating the matrix with transitions
  nstates <- length(states)
  mat <- matrix(
    0, nrow = nstates, ncol = nstates,
    dimnames = list(states, states)
  )

  for (s_i in states) {
    for (s_j in states) {
      # Filtering and aggregating
      mat[s_i, s_j] <- x[x$from == s_i & x$to == s_j, ]$counts |> sum()
    }
  }

  row_totals <- rowSums(mat)
  rows_with_transitions <- row_totals > 0
  if (any(rows_with_transitions)) {
    mat[rows_with_transitions, ] <-
      mat[rows_with_transitions, , drop = FALSE] /
        row_totals[rows_with_transitions]
  }

  draw_mermaid_from_matrix(mat, ...) |>
    plot()

}

#' @export
plot.epiworld_multiple_save_i <- function(x, y = NULL, ...) {

  what <- attr(x, "what")

  if (nrow(x) == 0) {
    warning(
      "When plotting the saved results, for the case of ",
      what, ", there were no observations."
    )
    return(NULL)
  }

  # If it is not reproductive number, then...
  if (what == "reproductive") {
    plot.epiworld_multiple_save_reproductive_number(x, ...)
  } else if (what == "transition") {
    plot_multiple_transition(x, ...)
  } else {

    states <- unique(x$state)
    ncols <- max(1L, ceiling(length(states) / 2))
    oldpar <- graphics::par(
      mfrow = c(2L, ncols)
    )
    on.exit(graphics::par(oldpar))

    for (state_name in states) {

      graphics::boxplot(
        counts ~ date,
        data = x[x$state == state_name, , drop = FALSE],
        main = state_name,
        xlab = "Date",
        ylab = "Counts",
        border = "black",
        las = 2
      )

    }

  }


}

#' @export
plot.epiworld_multiple_save_reproductive_number <- function(x, y = NULL, ...) {
  # Identifying sims
  sims <- sort(unique(x[["sim_num"]]))

  totals <- NULL
  for (s in sims) {
    # Subsetting the data
    x_tmp <- x[x[["sim_num"]] == s, , drop = FALSE]

    # Computing daily values
    totals <- rbind(
      totals,
      plot.epiworld_repnum(x_tmp, plot = FALSE)
    )

  }

  graphics::boxplot(
    avg ~ date,
    data = totals,
    main = "Reproductive Number",
    xlab = "Source Exposure Date",
    ylab = "rt",
    border = "black",
    las = 2
  )

  invisible(totals)

}
