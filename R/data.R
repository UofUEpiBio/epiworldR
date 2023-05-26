#' Accessing the database of `epiworld`
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
#'   prob_transmission = 0.8,
#'   incubation_days   = 7.0,
#'   prob_recovery     = 0.3
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
#' # Retrieving the reproductive number
#' get_reproductive_number(seirconn)
#' 
#' # We can go further and get all the history
#' t_hist <- get_hist_transition_matrix(seirconn)
#' 
#' head(t_hist)
#' 
#' # And turn it into an array
#' as.array(t_hist)[,,1:3]
NULL

#' @export
#' @rdname epiworld-data
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
#' @rdname epiworld-data
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
#' @rdname epiworld-data
get_reproductive_number <- function(x) UseMethod("get_reproductive_number")

#' @export
get_reproductive_number.epiworld_model <- function(x) {
  res <- get_reproductive_number_cpp(x)
  class(res) <- c("epiworld_repnum", class(res))
  res
}

#' @rdname epiworld-data
#' @param y Ignored.
#' @param plot Logical scalar. If `TRUE` (default), the function will plot the
#' average reproductive number over time.
#' @param ylab,xlab,main,type Further parameters passed to [graphics::plot()]
#' @export
plot.epiworld_repnum <- function(
    x,
    y = NULL,
    ylab = "Average Rep. Number",
    xlab = "Day (step)",
    main = "Daily Average Reproductive Number",
    plot = TRUE,
    ...) {
  
  
  # Nvariants
  vlabs     <- sort(unique(x[, "variant"]))
  nvariants <- length(vlabs)
  
  res <- vector("list", nvariants)
  names(res) <- vlabs
  for (i in seq_along(vlabs)) {
    
    x_tmp <- x[x[, "variant"] == vlabs[i], ]
    
    res[[i]] <- tapply(
      X     = x_tmp[, "rt"],
      INDEX = x_tmp[, "source_exposure_date"],
      FUN   = mean
    )
    
    # Preparing the data frame
    res[[i]] <- data.frame(
      variant = vlabs[i],
      step = as.integer(names(res[[i]])),
      avg  = unname(res[[i]])
    )
    
  }
  
  # Figuring out the range
  res_all <- do.call(rbind, res)
  
  yran <- range(res_all[["avg"]])
  xran <- range(res_all[["step"]])
  
  # Plotting -------------------------------------------------------------------
  if (plot) {
    for (i in seq_along(vlabs)) {
      
      if (i == 1L) {
        
        graphics::plot(
          x = res[[i]][["step"]],
          y = res[[i]][["avg"]],
          pch  = (i - 1),
          col  = i,
          xlab = xlab,
          ylab = ylab,
          main = main,
          ...
        )
        next
      }
      
      graphics::points(
        x = res[[i]][["step"]],
        y = res[[i]][["avg"]],
        pch  = (i - 1),
        col  = i,
        ...
      )
      
    }

    if (nvariants > 1L) {
      
      graphics::legend(
        "topright",
        legend = vlabs,
        pch    = 0L:(nvariants - 1L),
        col    = 1L:nvariants,
        title  = "Variants",
        bty    = "n"
      )
      
    }

  }
  
  res <- do.call(rbind, res)
  rownames(res) <- NULL
  invisible(res)
  
}

#' @export
#' @rdname epiworld-data
#' @returns 
#' In the case of `get_hist_transition_matrix`, a [data.frame] with four columns:
#' "state_from", "state_to", "date", and "counts." It will be
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
#' The `as.array` method for `epiworld_hist_transition` objects turns the `data.frame`
#' returned by `get_hist_transition_matrix` into an array of `nstates x nstates x (ndays + 1)`
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
#' @details The `plot_incidence` function is a wrapper between
#' [get_hist_transition_matrix] and it's plot method.
plot_incidence <- function(x, ...) {
  plot(get_hist_transition_matrix(x), ...)
}

#' @export 
#' @rdname epiworld-data
#' @details The plot method for the `epiworld_hist_transition` class plots the
#' daily incidence of each state. The function returns the data frame used for
#' plotting.
plot.epiworld_hist_transition <- function(
  x,
  type = "l",
  xlab = "Day (step)",
  ylab = "Counts",
  main = "Daily incidence",
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
  res <- res[, is_not_zero, drop = FALSE]
  states <- colnames(res)
  n_states <- length(states)

  # Plotting each column. We start by taking the total
  # range
  yran <- range(res)
  xran <- range(0:n_steps)
  for (i in seq_along(states)) {
    
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
  
  invisible(res)

}

#' @export
#' @rdname epiworld-data
#' @return The function `get_transmissions` returns a `data.frame` with the following
#' columns: `date`, `source`, `target`, `variant`, and `source_exposure_date`.
get_transmissions <- function(x) {
  
  stopifnot_model(x)
  res <- get_transmissions_cpp(x)
  structure(
    res,
    class = c("epiworld_transmissions", class(res))
  )
  
}