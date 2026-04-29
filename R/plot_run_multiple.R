#' Visualize results from `run_multiple()`
#' @param x An object of class `epiworld_run_multiple`
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

# Function to compute the 95% CI
calc_ci_ts <- function(
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

# ans2 <- calc_ci_ts(x = r$rt, group = r$source_exposure_date)
# ans3 <- calc_ci_ts(x = r$rt, group = r$source_exposure_date, alpha = .2)

# plot(
#   NULL,
#   xlim = range(ans2$area$x),
#   ylim = range(ans2$area$y)
# )
# polygon(ans2$area, col = "lightgrey", border = "darkgray")
# polygon(ans3$area, col = "darkgray", border = "black")
# lines(
#   ans2$line, type = "l",
#   lwd = 2, lty=2
#   )


#' @export
plot_epiworld_multiple_save_transition <- function(x, ...) {
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
      mat[s_i, s_j] <- subset(
        x, subset = from == s_i & to == s_j
      )$counts |> sum()
    }
  }

  mat <- mat / rowSums(mat)

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
    plot_epiworld_multiple_save_reproductive_number(x, ...)
  } else if (what == "transition") {
    plot_epiworld_multiple_save_transition(x, ...)
  } else {

    oldpar <- graphics::par(
      mfrow = c(2, floor(length(unique(x$state)) / 2))
    )
    on.exit(graphics::par(oldpar))

    for (what in unique(x$state)) {

      graphics::boxplot(
        counts ~ date,
        data = x[x$state == what, , drop = FALSE],
        main = what,
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
