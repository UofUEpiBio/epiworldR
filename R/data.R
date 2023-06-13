#' Accessing the database of `epiworld`
#' 
#' Models in `epiworld` are stored in a database. This database can be accessed
#' using the functions described in this manual page. Some elements of the
#' database are: the transition matrix, the incidence matrix, the reproductive 
#' number, the generation time, and daily incidence at the virus and tool level.
#' 
#' 
#' @param x An object of class [`epiworld_sir`], [`epiworld_seir`], etc.
#' any model.
#' @param skip_zeros Logical scalar. When `FALSE` it will return all the
#' entries in the transition matrix.
#' @param ... In the case of plot methods, further arguments passed to
#' [graphics::plot].
#' @name epiworld-data
#' @family Models
#' @examples
#' # SEIR Connected
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
#' # Retrieving the transition probability
#' get_transition_probability(seirconn) 
#' 
#' # Retrieving date, state, and counts dataframe including any added tools 
#' get_hist_tool(seirconn)
#' 
#' # Retrieving overall date, state, and counts dataframe 
#' head(get_hist_total(seirconn))
#' 
#' # Retrieving date, state, and counts dataframe by variant
#' head(get_hist_variant(seirconn))
#' 
#' # Retrieving (and plotting) the reproductive number
#' rp <- get_reproductive_number(seirconn)
#' plot(rp) # Also equivalent to plot_reproductive_number(seirconn)
#' 
#' # We can go further and get all the history
#' t_hist <- get_hist_transition_matrix(seirconn)
#' 
#' head(t_hist)
#' 
#' # And turn it into an array
#' as.array(t_hist)[,,1:3]
#' 
#' # We cam also get (and plot) the incidence, as well as
#' # the generation time
#' inci <- plot_incidence(seirconn)
#' gent <- plot_generation_time(seirconn)
#' 
NULL

#' @export
#' @returns 
#' - The `get_hist_total` function returns an object of class 
#' [epiworld_hist_total].
#' @rdname epiworld-data
#' @aliases epiworld_hist_total
get_hist_total <- function(x) UseMethod("get_hist_total")

#' @export
get_hist_total.epiworld_model <- function(x)  {
  
  res <- get_hist_total_cpp(x)
  structure(
    res,
    class  = c("epiworld_hist_total", "epiworld_hist", "data.frame"),
    states = sort(unique(res$state))
  )
  
}

#' @export
plot.epiworld_hist <- function(x, y, ...) {
  plot_epi(x, ...)
}

#' @export
#' @returns 
#' - The `get_hist_variant` function returns an object of class 
#' [epiworld_hist_variant].
#' @rdname epiworld-data
#' @aliases epiworld_hist_variant
get_hist_variant <- function(x) UseMethod("get_hist_variant")

#' @export
get_hist_variant.epiworld_model <- function(x)  {
  res <- get_hist_variant_cpp(x)
  
  structure(
    res,
    class  = c("epiworld_hist_variant", "epiworld_hist", "data.frame"),
    states = sort(unique(res$state))
  )
}

#' @export
#' @returns 
#' - The `get_hist_tool` function returns an object of [epiworld_hist_variant].
#' @rdname epiworld-data
get_hist_tool <- function(x) UseMethod("get_hist_tool")

#' @export
get_hist_tool.epiworld_model <- function(x)  {
  res <- get_hist_tool_cpp(x)
  structure(
    res,
    class  = c("epiworld_hist_tool", "epiworld_hist", "data.frame"),
    states = sort(unique(res$state))
  )
}

#' @export
#' @returns 
#' - The `get_transition_probability` function returns an object of class 
#' `matrix`. 
#' @rdname epiworld-data
get_transition_probability <- function(x) {
  UseMethod("get_transition_probability")
}

#' @export
get_transition_probability.epiworld_model <- function(x)  {
  res <- get_transition_probability_cpp(x)
  s   <- get_states(x)
  
  ns <- length(s)
  
  matrix(res, nrow = ns, ncol = ns, dimnames = list(s, s))
}

#' @export
#' @returns 
#' - The `get_reproductive_number` function returns an object of class
#' [epiworld_repnum].
#' @rdname epiworld-data
#' @aliases epiworld_repnum
get_reproductive_number <- function(x) UseMethod("get_reproductive_number")

#' @export
get_reproductive_number.epiworld_model <- function(x) {
  res <- get_reproductive_number_cpp(x)
  class(res) <- c("epiworld_repnum", class(res))
  res
}

#' @rdname epiworld-data
#' @param y Ignored.
#' @param plot Logical scalar. If `TRUE` (default), the function will the
#' desired statistic.
#' @param ylab,xlab,main,type Further parameters passed to [graphics::plot()]
#' @returns 
#' - The `plot` function returns a plot of the reproductive number over time. 
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
      variant = integer(),
      date    = integer(),
      avg     = numeric(),
      n       = integer(),
      sd      = numeric(),
      lb      = numeric(),
      ub      = numeric()
    )
  } else {

    # Computing stats
    # Compute the mean and 95% CI of rt by variant and source_exposure_date using the repnum data.frame with the tapply function
    repnum <- stats::aggregate(
      x[["rt"]],
      by = list(
        variant = x[["variant"]],
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
    repnum <- repnum[order(repnum[["variant"]], repnum[["date"]]), , drop = FALSE]
    rownames(repnum) <- NULL
  }

  # Nvariants
  vlabs     <- sort(unique(x[, "variant"]))
  nvariants <- length(vlabs)
    
  # # Figuring out the range
  yran <- range(repnum[["avg"]], na.rm = TRUE)
  xran <- range(repnum[["date"]], na.rm = TRUE)
  
  # Plotting -------------------------------------------------------------------
  if (plot) {
    for (i in seq_along(vlabs)) {

      tmp <- repnum[repnum[["variant"]] == vlabs[i], ]
      
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

    if (nvariants > 1L) {
      
      graphics::legend(
        "topright",
        legend = vlabs,
        pch    = 1L:nvariants,
        col    = 1L:nvariants,
        lwd    = 2,
        title  = "Variants",
        bty    = "n"
      )
      
    }

  }
  invisible(repnum)
  
}

#' @export
#' @rdname epiworld-data
#' @details The `plot_reproductive_number` function is a wrapper around
#' [get_reproductive_number] that plots the result.
plot_reproductive_number <- function(x, ...) {
  plot(get_reproductive_number(x), ...)
}



#' @export
#' @rdname epiworld-data
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
#' @rdname epiworld-data
as.array.epiworld_hist_transition <- function(x, ...) {
  
  states <- attr(x, "states")
  n_states <- length(states)
  n_steps  <- attr(x, "nsteps")
  
  res <- array(
    0L,
    dim      = c(n_states, n_states, n_steps + 1), # Includes the baseline
    dimnames = list(states, states, 0:n_steps)
    )
 
  res[cbind(x[,1], x[,2], x[,3])] <- x[,4]
  
  res

}

#' @export
#' @rdname epiworld-data
#' @returns 
#' - The `plot_incidence` function returns a plot originating from the object 
#' `get_hist_transition_matrix`.
#' @details The `plot_incidence` function is a wrapper between
#' [get_hist_transition_matrix] and it's plot method.
plot_incidence <- function(x, ...) {
  plot(get_hist_transition_matrix(x), ...)
}

#' @export 
#' @returns
#' - The `plot` function returns a plot which originates from the 
#' `epiworld_hist_transition` object.
#' @rdname epiworld-data
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

#' @export
#' @rdname epiworld-data
#' @return
#' - The function `get_transmissions` returns a `data.frame` with the following
#' columns: `date`, `source`, `target`, `variant`, and `source_exposure_date`.
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

#' @export
#' @rdname epiworld-data
#' @return
#' The function `get_generation_time` returns a `data.frame` with
#' the following columns: "agent", "virus_id", "date", and "gentime".
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
#' @rdname epiworld-data
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

  gt <- stats::aggregate(
    x[["gentime"]], by = list(
      date    = x[["date"]],
      variant = x[["virus_id"]]
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
  gt <- gt[order(gt[["variant"]], gt[["date"]]), , drop = FALSE]
  rownames(gt) <- NULL

  # Replacing NaNs with NAs
  gt <- as.data.frame(lapply(gt, function(x) {
    x[is.nan(x)] <- NA
    x
  })) 

  if (plot) {
    # Number of viruses
    variants <- sort(unique(gt[["variant"]]))
    n_viruses <- length(variants)

    for (i in 1L:n_viruses) {

      gt_i <- gt[gt[["variant"]] == variants[i], , drop = FALSE]
      
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
        legend = variants,
        col    = 1:n_viruses,
        lwd    = 2,
        lty    = 1L:n_viruses,
        title  = "Viruses",
        bty    = "n"
      )
      
    }

  }

  invisible(gt)
    
  
}

#' @export 
#' @rdname epiworld-data
#' @return 
#' - The function `plot_generation_time` is a wrapper for [plot] and
#' [get_generation_time].
plot_generation_time <- function(x, ...) {
  plot(get_generation_time(x), ...)
}

