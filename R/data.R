# ============================================================================
# TOPIC: epiworld-history
# Model history and totals
# ============================================================================

#' Model history and totals
#'
#' Functions to extract simulation history at total, variant, and tool levels,
#' plus snapshot totals and a common plot method for history objects.
#'
#' @param x An object of class [`epiworld_sir`], [`epiworld_seir`], etc.
#' (any model).
#' @param ... In the case of plot methods, further arguments passed to
#' [graphics::plot].
#' @name epiworld-history
#' @concept model-utility-functions
#' @family History
#' @examples
#' # SEIR Connected model
#' seirconn <- ModelSEIRCONN(
#'   name              = "Disease",
#'   n                 = 10000,
#'   prevalence        = 0.1,
#'   contact_rate      = 2.0,
#'   transmission_rate = 0.8,
#'   incubation_days   = 7.0,
#'   recovery_rate     = 0.3
#' )
#'
#' # Running the simulation for 50 steps (days)
#' set.seed(937)
#' run(seirconn, 50)
#'
#' # Retrieving date, state, and counts dataframe including any added tools
#' get_hist_tool(seirconn)
#'
#' # Retrieving overall date, state, and counts dataframe
#' head(get_hist_total(seirconn))
#'
#' # Retrieving date, state, and counts dataframe by variant
#' head(get_hist_virus(seirconn))
#'
#' # Snapshot of totals at end of simulation
#' get_today_total(seirconn)
#'
NULL

#' @export
#' @returns
#' - The `get_hist_total` function returns an object of class
#' [epiworld_hist_total].
#' @rdname epiworld-history
#' @aliases epiworld_hist_total
get_hist_total <- function(x) UseMethod("get_hist_total")

#' @export
get_hist_total.epiworld_model <- function(x) {

  res <- get_hist_total_cpp(x)
  structure(
    res,
    class  = c("epiworld_hist_total", "epiworld_hist", "data.frame"),
    states = sort(unique(res$state))
  )

}

#' @export
#' @rdname epiworld-history
#' @returns
#' - The `get_today_total` function returns a named vector with the
#' total number of individuals in each state at the end of the simulation.
get_today_total <- function(x) UseMethod("get_today_total")

#' @export
get_today_total.epiworld_model <- function(x) {
  get_today_total_cpp(x)
}

#' @export
#' @rdname epiworld-history
#' @param y Ignored.
plot.epiworld_hist <- function(x, y, ...) {
  plot_epi(x, ...)
}

#' @export
#' @returns
#' - The `get_hist_virus` function returns an object of class
#' [epiworld_hist_virus].
#' @rdname epiworld-history
#' @aliases epiworld_hist_variant epiworld_hist_virus
get_hist_virus <- function(x) UseMethod("get_hist_virus")

#' @export
get_hist_virus.epiworld_model <- function(x) {
  res <- get_hist_virus_cpp(x)

  structure(
    res,
    class  = c("epiworld_hist_virus", "epiworld_hist", "data.frame"),
    states = sort(unique(res$state))
  )
}

#' @export
#' @returns
#' - The `get_hist_tool` function returns an object of [epiworld_hist_tool].
#' @rdname epiworld-history
#' @aliases epiworld_hist_tool
get_hist_tool <- function(x) UseMethod("get_hist_tool")

#' @export
get_hist_tool.epiworld_model <- function(x) {
  res <- get_hist_tool_cpp(x)
  structure(
    res,
    class  = c("epiworld_hist_tool", "epiworld_hist", "data.frame"),
    states = sort(unique(res$state))
  )
}

# ============================================================================
# TOPIC: epiworld-summaries
# Summary counts and probabilities
# ============================================================================

#' Summary counts and probabilities
#'
#' Functions to extract summary statistics from models, including transition
#' probabilities, active cases, and outbreak sizes.
#'
#' @param x An object of class [`epiworld_sir`], [`epiworld_seir`], etc.
#' (any model).
#' @name epiworld-summaries
#' @concept model-utility-functions
#' @family Summaries
#' @examples
#' # SEIR Connected model
#' seirconn <- ModelSEIRCONN(
#'   name              = "Disease",
#'   n                 = 10000,
#'   prevalence        = 0.1,
#'   contact_rate      = 2.0,
#'   transmission_rate = 0.8,
#'   incubation_days   = 7.0,
#'   recovery_rate     = 0.3
#' )
#'
#' set.seed(937)
#' run(seirconn, 50)
#'
#' # Retrieving the transition probability
#' get_transition_probability(seirconn)
#'
#' # Get active cases
#' head(get_active_cases(seirconn))
#'
#' # Get outbreak size
#' head(get_outbreak_size(seirconn))
#'
NULL

#' @export
#' @returns
#' - The `get_transition_probability` function returns an object of class
#' `matrix`.
#' @rdname epiworld-summaries
get_transition_probability <- function(x) {
  UseMethod("get_transition_probability")
}

#' @export
get_transition_probability.epiworld_model <- function(x) {
  res <- get_transition_probability_cpp(x)
  s   <- get_states(x)

  ns <- length(s)

  matrix(res, nrow = ns, ncol = ns, dimnames = list(s, s))
}

# ============================================================================
# TOPIC: epiworld-repnum
# Reproductive number (Rt)
# ============================================================================

#' Reproductive number (Rt)
#'
#' Extraction and plotting of the reproductive number (Rt) by virus over time.
#'
#' @param x An object of class [`epiworld_sir`], [`epiworld_seir`], etc.
#' (any model), or an object of class [epiworld_repnum].
#' @param ... In the case of plot methods, further arguments passed to
#' [graphics::plot].
#' @name epiworld-repnum
#' @concept model-utility-functions
#' @family Epidemiological metrics
#' @examples
#' # SEIR Connected model
#' seirconn <- ModelSEIRCONN(
#'   name              = "Disease",
#'   n                 = 10000,
#'   prevalence        = 0.1,
#'   contact_rate      = 2.0,
#'   transmission_rate = 0.8,
#'   incubation_days   = 7.0,
#'   recovery_rate     = 0.3
#' )
#'
#' set.seed(937)
#' run(seirconn, 50)
#'
#' # Retrieving (and plotting) the reproductive number
#' rp <- get_reproductive_number(seirconn)
#' plot(rp) # Also equivalent to plot_reproductive_number(seirconn)
#'
NULL

#' @export
#' @returns
#' - The `get_reproductive_number` function returns an object of class
#' [epiworld_repnum].
#' @rdname epiworld-repnum
#' @aliases epiworld_repnum
get_reproductive_number <- function(x) UseMethod("get_reproductive_number")

#' @export
get_reproductive_number.epiworld_model <- function(x) {
  res <- get_reproductive_number_cpp(x)
  res <- res[res[["source"]] != -1, , drop = FALSE]
  class(res) <- c("epiworld_repnum", class(res))
  res
}

#' @rdname epiworld-repnum
#' @param y Ignored.
#' @param plot Logical scalar. If `TRUE` (default), the function will plot the
#' desired statistic.
#' @param ylab,xlab,main,type Further parameters passed to [graphics::plot()]
#' @returns
#' - The `plot` method for `epiworld_repnum` returns a plot of the reproductive
#' number over time.
#' @export
#' @importFrom stats sd quantile aggregate
plot.epiworld_repnum <- function(
  x,
  y = NULL,
  ylab = "Average Rep. Number",
  xlab = "Day (step)",
  main = "Reproductive Number",
  type = "b",
  plot = TRUE,
  ...) {


  if (nrow(x) == 0) {
    repnum <- data.frame(
      virus_id = integer(),
      virus    = character(),
      date     = integer(),
      avg      = numeric(),
      n        = integer(),
      sd       = numeric(),
      lb       = numeric(),
      ub       = numeric()
    )
  } else {
    # Computing stats
    # Compute the mean and 95% CI of rt by virus and source_exposure_date using the repnum data.frame with the tapply function

    # Creating a new column combining virus_id and variant
    x[["virus_comb"]] <- sprintf("%s (%i)", x[["virus"]], x[["virus_id"]])

    repnum <- stats::aggregate(
      x[["rt"]],
      by = list(
        virus_comb = x[["virus_comb"]],
        date    = x[["source_exposure_date"]]
      ),
      FUN = function(x) {
        ci <- stats::quantile(x, c(0.025, 0.975), na.rm = TRUE)

        data.frame(
          avg  = mean(x, na.rm = TRUE),
          n    = sum(!is.na(x)),
          sd   = stats::sd(x, na.rm = TRUE),
          lb   = ci[1],
          ub   = ci[2]
        )
      },
      simplify = FALSE
    )

    repnum <- cbind(repnum[, -3, drop = FALSE], do.call(rbind, repnum[, 3]))
    repnum <- repnum[order(repnum[["virus_comb"]], repnum[["date"]]), , drop = FALSE]

    # Merging the virus and virus_id column of x to repnum
    repnum <- merge(
      repnum,
      unique(x[, c("virus", "virus_id", "virus_comb")]),
      by = "virus_comb",
      all.x = TRUE,
      all.y = FALSE
    )

    rownames(repnum) <- NULL

    # Reordering columns
    repnum <- repnum[, c(
      "virus_id", "virus", "date", "avg", "n", "sd", "lb", "ub",
      "virus_comb"
    )]

  }


  # Nviruses
  vlabs     <- sort(unique(x[, "virus_comb"]))
  nviruses <- length(vlabs)

  # # Figuring out the range
  yran <- range(repnum[["avg"]], na.rm = TRUE)
  xran <- range(repnum[["date"]], na.rm = TRUE)

  # Plotting -------------------------------------------------------------------
  if (plot) {
    for (i in seq_along(vlabs)) {

      tmp <- repnum[repnum[["virus_comb"]] == vlabs[i], ]

      if (i == 1L) {

        graphics::plot(
          x = tmp[["date"]],
          y = tmp[["avg"]],
          pch  = i,
          col  = i,
          lwd  = 2,
          lty  = i,
          xlab = xlab,
          ylab = ylab,
          main = main,
          type = type,
          xlim = xran,
          ylim = yran,
          ...
        )
        next
      }

      graphics::lines(
        x = tmp[["date"]],
        y = tmp[["avg"]],
        pch  = i,
        col  = i,
        lwd  = 2,
        lty  = i,
        type = type,
        ...
      )

    }

    if (nviruses > 1L) {

      graphics::legend(
        "topright",
        legend = vlabs,
        pch    = 1L:nviruses,
        col    = 1L:nviruses,
        lwd    = 2,
        title  = "Virus",
        bty    = "n"
      )

    }

  }

  # Removing the virus_comb column
  repnum[["virus_comb"]] <- NULL

  invisible(repnum)

}

#' @export
#' @rdname epiworld-repnum
#' @details The `plot_reproductive_number` function is a wrapper around
#' [get_reproductive_number] that plots the result.
plot_reproductive_number <- function(x, ...) {
  plot(get_reproductive_number(x), ...)
}


# ============================================================================
# TOPIC: epiworld-transition
# Transition dynamics and incidence
# ============================================================================

#' Transition dynamics and incidence
#'
#' Functions to extract and visualize state transition counts, daily incidence,
#' and conversion to array format.
#'
#' @param x An object of class [`epiworld_sir`], [`epiworld_seir`], etc.
#' (any model), or an object of class [epiworld_hist_transition].
#' @param skip_zeros Logical scalar. When `FALSE` it will return all the
#' entries in the transition matrix.
#' @param ... In the case of plot methods, further arguments passed to
#' [graphics::plot].
#' @name epiworld-transition
#' @aliases epiworld_hist_transition
#' @concept model-utility-functions
#' @family Transition dynamics
#' @examples
#' # SEIR Connected model
#' seirconn <- ModelSEIRCONN(
#'   name              = "Disease",
#'   n                 = 10000,
#'   prevalence        = 0.1,
#'   contact_rate      = 2.0,
#'   transmission_rate = 0.8,
#'   incubation_days   = 7.0,
#'   recovery_rate     = 0.3
#' )
#'
#' set.seed(937)
#' run(seirconn, 50)
#'
#' # Get the transition history
#' t_hist <- get_hist_transition_matrix(seirconn)
#' head(t_hist)
#'
#' # Convert to array
#' as.array(t_hist)[, , 1:3]
#'
#' # Plot incidence
#' inci <- plot_incidence(seirconn)
#'
NULL

#' @export
#' @rdname epiworld-transition
#' @returns
#' - `get_hist_transition_matrix` returns a [data.frame] with four columns:
#' "state_from", "state_to", "date", and "counts."
get_hist_transition_matrix <- function(x, skip_zeros = FALSE)
  UseMethod("get_hist_transition_matrix")

#' @export
get_hist_transition_matrix.epiworld_model <- function(x, skip_zeros = FALSE) {

  res <- get_hist_transition_matrix_cpp(x, skip_zeros)
  class(res) <- c("epiworld_hist_transition", class(res))

  attr(res, "states") <- get_states(x)
  attr(res, "nsteps") <- get_ndays(x)

  res

}


#' @export
#' @returns
#' - The `as.array` method for `epiworld_hist_transition` objects turns the
#' `data.frame` returned by `get_hist_transition_matrix` into an array of
#' `nstates x nstates x (ndays + 1)`
#' entries, where the first entry is the initial state.
#' @rdname epiworld-transition
as.array.epiworld_hist_transition <- function(x, ...) {

  states <- attr(x, "states")
  n_states <- length(states)
  n_steps  <- attr(x, "nsteps")

  res <- array(
    0L,
    dim      = c(n_states, n_states, n_steps + 1), # Includes the baseline
    dimnames = list(states, states, 0:n_steps)
  )

  res[cbind(x[, 1], x[, 2], x[, 3])] <- x[, 4]

  res

}

#' @export
#' @rdname epiworld-transition
#' @returns
#' - The `plot_incidence` function returns a plot originating from the object
#' `get_hist_transition_matrix`.
#' @details The `plot_incidence` function is a wrapper between
#' [get_hist_transition_matrix] and its plot method.
plot_incidence <- function(x, ...) {
  plot(get_hist_transition_matrix(x), ...)
}

#' @export
#' @returns
#' - The `plot` method for `epiworld_hist_transition` returns a plot of the
#' daily incidence.
#' @rdname epiworld-transition
#' @param plot Logical scalar. If `TRUE` (default), the function will plot the
#' desired statistic.
#' @param ylab,xlab,main,type Further parameters passed to [graphics::plot()]
#' @details The plot method for the `epiworld_hist_transition` class plots the
#' daily incidence of each state. The function returns the data frame used for
#' plotting.
plot.epiworld_hist_transition <- function(
  x,
  type = "b",
  xlab = "Day (step)",
  ylab = "Counts",
  main = "Daily incidence",
  plot = TRUE,
  ...
) {

  if (!inherits(x, "epiworld_hist_transition")) {
    stop("The object must be of class 'epiworld_hist_transition'")
  }

  states <- attr(x, "states")
  n_states <- length(states)
  n_steps  <- attr(x, "nsteps")

  # Agregating the data
  x <- x[x[, "state_from"] != x[, "state_to"], , drop = FALSE]
  res <- tapply(x[["counts"]], INDEX = list(x[, "state_to"], x[, "date"]), FUN = sum)
  res <- as.data.frame(t(res))

  # Checking if any of the columns is all zeros
  is_not_zero <- which(colSums(res) != 0)
  # res <- res[, is_not_zero, drop = FALSE]
  states <- colnames(res)
  n_states <- length(states)

  # Plotting each column. We start by taking the total
  # range
  if (plot) {
    yran <- range(res)
    xran <- range(0:n_steps)
    for (i in is_not_zero) {

      col <- states[i]

      if (i == 1L) {

        graphics::plot(
          x = as.integer(rownames(res)),
          y = res[[col]],
          col = i,
          lwd = 2,
          lty = i,
          type = type,
          xlab = xlab,
          ylab = ylab,
          main = main,
          ylim = yran,
          xlim = xran,
          ...
        )

        next

      }

      graphics::points(
        x = as.integer(rownames(res)),
        y = res[[col]],
        type = type,
        col = i,
        lwd = 2,
        lty = i,
        ...
      )

    }

    # Creating a legend
    if (n_states > 1L) {

      graphics::legend(
        "topright",
        legend = states,
        col    = 1L:n_states,
        lwd    = 2,
        lty    = 1L:n_states,
        title  = "States",
        bty    = "n"
      )

    }
  }

  invisible(res)

}

# ============================================================================
# TOPIC: epiworld-transmissions
# Transmission network
# ============================================================================

#' Transmission network
#'
#' Transmission edges, including seeded infections (`source = -1`).
#'
#' @param x An object of class [`epiworld_sir`], [`epiworld_seir`], etc.
#' (any model).
#' @name epiworld-transmissions
#' @concept model-utility-functions
#' @family Network outputs
#' @examples
#' # SEIR Connected model
#' seirconn <- ModelSEIRCONN(
#'   name              = "Disease",
#'   n                 = 10000,
#'   prevalence        = 0.1,
#'   contact_rate      = 2.0,
#'   transmission_rate = 0.8,
#'   incubation_days   = 7.0,
#'   recovery_rate     = 0.3
#' )
#'
#' set.seed(937)
#' run(seirconn, 50)
#'
#' # Get transmission data
#' head(get_transmissions(seirconn))
#'
NULL

#' @export
#' @rdname epiworld-transmissions
#' @details
#' The function `get_transmissions` includes the seeded infections, with the
#' `source` column coded as `-1`.
#' @return
#' - The function `get_transmissions` returns a `data.frame` with the following
#' columns: `date`, `source`, `target`, `virus_id`, `virus`, and `source_exposure_date`.
get_transmissions <- function(x) UseMethod("get_transmissions")

#' @export
get_transmissions.epiworld_diffnet <- function(x) {
  warning("The transmission network is not necesarily relevant for the diffnet model")
  get_transmissions.epiworld_model(x)
}

#' @export
get_transmissions.epiworld_model <- function(x) {

  res <- get_transmissions_cpp(x)
  structure(
    res,
    class = c("epiworld_transmissions", class(res))
  )

}

# ============================================================================
# TOPIC: epiworld-gentime
# Generation time
# ============================================================================

#' Generation time
#'
#' Extraction and plotting of generation time by virus over time.
#'
#' @param x An object of class [`epiworld_sir`], [`epiworld_seir`], etc.
#' (any model), or an object of class [epiworld_generation_time].
#' @param ... In the case of plot methods, further arguments passed to
#' [graphics::plot].
#' @name epiworld-gentime
#' @aliases epiworld_generation_time
#' @concept model-utility-functions
#' @family Epidemiological metrics
#' @examples
#' # SEIR Connected model
#' seirconn <- ModelSEIRCONN(
#'   name              = "Disease",
#'   n                 = 10000,
#'   prevalence        = 0.1,
#'   contact_rate      = 2.0,
#'   transmission_rate = 0.8,
#'   incubation_days   = 7.0,
#'   recovery_rate     = 0.3
#' )
#'
#' set.seed(937)
#' run(seirconn, 50)
#'
#' # Get and plot generation time
#' gent <- plot_generation_time(seirconn)
#'
NULL

#' @export
#' @rdname epiworld-gentime
#' @return
#' - The function `get_generation_time` returns a `data.frame` with
#' the following columns: "agent", "virus_id", "virus", "date", and "gentime".
get_generation_time <- function(x) {

  stopifnot_model(x)
  res <- get_generation_time_cpp(x)

  # Replacing -1 with NAs
  res[["gentime"]][res[["gentime"]] == -1] <- NA_integer_

  structure(
    res,
    class = c("epiworld_generation_time", class(res)),
    n_steps = get_ndays(x)
  )

}

#' @export
#' @rdname epiworld-gentime
#' @param plot Logical scalar. If `TRUE` (default), the function will plot the
#' desired statistic.
#' @param ylab,xlab,main,type Further parameters passed to [graphics::plot()]
plot.epiworld_generation_time <- function(
  x,
  type = "b",
  xlab = "Day (step)",
  ylab = "Avg. Generation Time",
  main = "Generation Time",
  plot = TRUE,
  ...
) {

  if (!inherits(x, "epiworld_generation_time")) {
    stop("The object must be of class 'epiworld_generation_time'")
  }

  # Combining virus with virus id (as we've done before)
  x[["virus_comb"]] <- sprintf(
    "%s (%s)",
    x[["virus"]],
    x[["virus_id"]]
  )

  gt <- stats::aggregate(
    x[["gentime"]], by = list(
      date    = x[["date"]],
      virus_comb = x[["virus_comb"]]
    ),
    FUN = function(x) {
      ci <- stats::quantile(
        x, probs = c(0.025, 0.975), na.rm = TRUE
      )

      data.frame(
        avg = mean(x, na.rm = TRUE),
        n   = sum(!is.na(x)),
        sd  = sd(x, na.rm = TRUE),
        ci_lower = ci[1L],
        ci_upper = ci[2L]
      )
    },
    simplify = FALSE
  )

  gt <- cbind(gt[, -3, drop = FALSE], do.call(rbind, gt[, 3]))
  gt <- gt[order(gt[["virus_comb"]], gt[["date"]]), , drop = FALSE]

  # Merging the virus and virus_id column of x to repnum
  gt <- merge(
    gt,
    unique(x[, c("virus", "virus_id", "virus_comb")]),
    by = "virus_comb",
    all.x = TRUE,
    all.y = FALSE
  )

  rownames(gt) <- NULL

  # Replacing NaNs with NAs
  gt <- as.data.frame(lapply(gt, function(x) {
    x[is.nan(x)] <- NA
    x
  }))

  if (plot) {
    # Number of viruses
    viruses <- sort(unique(gt[["virus_comb"]]))
    n_viruses <- length(viruses)

    for (i in 1L:n_viruses) {

      gt_i <- gt[gt[["virus_comb"]] == viruses[i], , drop = FALSE]

      if (i == 1L) {

        graphics::plot(
          x = gt_i[["date"]],
          y = gt_i[["avg"]],
          col = i,
          lwd = 2,
          lty = i,
          type = type,
          xlab = xlab,
          ylab = ylab,
          main = main,
          ylim = range(gt[["avg"]], na.rm = TRUE),
          xlim = range(gt[["date"]], na.rm = TRUE),
          ...
        )

        next

      }

      graphics::points(
        x = gt_i[["date"]],
        y = gt_i[["avg"]],
        col = i,
        lwd = 2,
        lty = i,
        type = type,
        ...
      )
    }

    # Creating a legend
    if (n_viruses > 1L) {

      graphics::legend(
        "topright",
        legend = viruses,
        col    = 1:n_viruses,
        lwd    = 2,
        lty    = 1L:n_viruses,
        title  = "Virus",
        bty    = "n"
      )

    }

  }

  # Deleting the virus_comb column
  gt[["virus_comb"]] <- NULL

  invisible(gt)


}

#' @export
#' @rdname epiworld-gentime
#' @return
#' - The function `plot_generation_time` is a wrapper for [plot] and
#' [get_generation_time].
plot_generation_time <- function(x, ...) {
  plot(get_generation_time(x), ...)
}

#' @export
#' @rdname epiworld-summaries
#' @return
#' - The function `get_active_cases` returns a data.frame with four columns:
#' date, virus_id, virus, and active_cases indicating the number of active
#' cases (individuals with a virus) at each point in time.
get_active_cases <- function(x) {
  stopifnot_model(x)
  return(get_active_cases_cpp(x))
}

#' @export
#' @rdname epiworld-summaries
#' @return
#' - The function `get_outbreak_size` returns a data.frame with four columns:
#' `date`, `virus_id`, `virus`, and `outbreak_size` indicating the outbreak
#' size per virus at each point in time.
get_outbreak_size <- function(x) {
  stopifnot_model(x)
  return(get_outbreak_size_cpp(x))
}

# ============================================================================
# TOPIC: epiworld-hospitalizations
# Hospitalizations by tool
# ============================================================================

#' Hospitalizations by tool
#'
#' Weighted hospitalization tracking when agents have multiple tools.
#'
#' @param x An object of class [`epiworld_sir`], [`epiworld_seir`], etc.
#' (any model).
#' @name epiworld-hospitalizations
#' @concept model-utility-functions
#' @family Summaries
#' @examples
#' # See model documentation for examples with hospitalization tracking
NULL

#' @export
#' @rdname epiworld-hospitalizations
#' @return
#' - The function `get_hospitalizations` returns a data.frame with five columns:
#' `date`, `virus_id`, `tool_id`, `hospitalizations`, and `weight`. The `weight`
#' column is used to keep track of individuals having multiple tools. For
#' example, if an agent has two tools (vaccination and mask-wearing), then it
#' will show up twice under count, but with weights 0.5 for each count. Models
#' with no hospitalization tracking will return the same data.frame with no
#' rows.
get_hospitalizations <- function(x) {

  stopifnot_model(x)
  return(get_hospitalizations_cpp(x))

}
