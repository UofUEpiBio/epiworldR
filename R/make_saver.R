#' Run multiple simulations at once
#' 
#' The `run_multiple` function allows running multiple simulations at once.
#' When available, users can take advantage of parallel computing to speed up
#' the process.
#' 
#' @param m,ndays,seed See [run].
#' @param saver An object of class [epiworld_saver].
#' @param nsims Integer. Number of replicats
#' @param ... List of strings (characters) specifying what to save (see details).
#' @param fn A file name pattern.
#' @param nthreads Integer. Number of threads (parallel computing.)
#' @param reset When `TRUE` (default,) resets the simulation.
#' @param verbose When `TRUE` (default,) prints a progress bar.
#' 
#' @details
#' Currently, the following elements can be saved:
#' 
#' - `total_hist` History of the model (total numbers per time).
#' - `variant_info` Information about `variants`.
#' - `variant_hist` Changes in `variants`.
#' - `tool_info` Information about `tools`.
#' - `tool_hist` Changes in `tools`.
#' - `transmission` Transmission events.
#' - `transition` Transition matrices.
#' - `reproductive` Reproductive number.
#' - `generation` Estimation of generation time.
#' 
#' @returns
#' - In the case of `make_saver`, an list of class `epiworld_saver`.
#' @examples
#' model_sir <- ModelSIRCONN(
#'   name = "COVID-19",
#'   prevalence = 0.01,
#'   n = 1000,
#'   contact_rate = 2,
#'   prob_transmission = 0.9, prob_recovery = 0.1
#'   )
#' 
#' # Generating a saver
#' saver <- make_saver("total_hist", "reproductive")
#' 
#' # Running and printing
#' run_multiple(model_sir, ndays = 100, nsims = 50, saver = saver, nthread = 2)
#' 
#' # Retrieving the results
#' ans <- run_multiple_get_results(model_sir)
#' 
#' head(ans$total_hist)
#' head(ans$reproductive)
#' 
#' # Plotting
#' multi_sir <- run_multiple_get_results(model_sir)$total_hist
#' multi_sir <- multi_sir[multi_sir$date <= 20,]
#' plot(multi_sir)
#' 
#' @export
#' @returns 
#' - The `run_multiple` function runs a specified number of simulations and 
#' returns a model object of class [epiworld_model].
run_multiple <- function(
    m, ndays, nsims,
    seed = sample.int(1e4, 1),
    saver = make_saver(),
    reset = TRUE,
    verbose = TRUE,
    nthreads = 1L
) UseMethod("run_multiple")

#' @export
run_multiple.epiworld_model <- function(
    m, ndays, nsims,
    seed     = sample.int(1e4, 1),
    saver    = make_saver(),
    reset    = TRUE,
    verbose  = TRUE,
    nthreads = 1L
) {
  
  if (!inherits(saver, "epiworld_saver"))
    stop("-saver- should be of class \"epiworld_saver\"")
  
  run_multiple_cpp(
    m,
    ndays,
    nsims,
    seed,
    saver$ptr,
    reset,
    verbose,
    nthreads
  )
  
  attr(m, "saver") <- saver
  
  invisible(m)
  
}

#' @export
#' @rdname run_multiple
#' @returns 
#' - The `run_multiple_get_results` function returns a named list with the
#' data specified by `make_saver`.
#' @importFrom utils read.table
run_multiple_get_results <- function(m) {
  
  if (!inherits(m, "epiworld_model"))
    stop("-m- must be of class `epiworld_model`.")
  
  # Get the filepath
  saver <- attr(m, "saver")
  
  if (!length(saver)) 
    stop("No -saver- found. -run_multiple_get_results- can only be used after using -run_multiple-.")
  
  pattern <- gsub("%[0-9]*lu", "*", saver$fn)
  
  output <- vector("list", length(saver$what))
  names(output) <- saver$what
  
  for (i in saver$what) {
    
    # Listing the files
    fnames <- list.files(
      path    = dirname(pattern),
      pattern = sprintf("%s\\.csv", i),
      full.names = TRUE
    )
    
    # Reading the files
    output[[i]] <- lapply(fnames, utils::read.table, sep = " ", header = TRUE)
    
    # Getting number of simulation
    output[[i]] <- lapply(seq_along(fnames), function(j) {
      if (nrow(output[[i]][[j]]) > 0)
        cbind(sim_num = j, output[[i]][[j]])
      else
        NULL
    })
    
    # Putting all together
    output[[i]] <- do.call(rbind, output[[i]])
    
    # runif
    # 
    # runif(10, min=2, max=4)
    # do.call(runif, list(10, 2, 4))
    # 
    # rbind(output[[i]][[1]], output[[i]][[2]], ...)
    
    class(output[[i]]) <- c("epiworld_multiple_save_i", class(output[[i]]))
    attr(output[[i]], "what") <- i
      
    
  }
  
  structure(output, class = c("epiworld_multiple_save", class(output)))
  
}

#' @export
plot.epiworld_multiple_save <- function(x, y = NULL, ...) {

  # what <- attr(x, "what")
  lapply(x, plot)

}

#' @export
plot.epiworld_multiple_save_i <- function(x, y = NULL, ...) {

  what <- attr(x, "what")
  
  # If it is not reproductive number, then...
  if (what != "reproductive") {
    
    oldpar <- graphics::par(
      mfrow = c(2, floor(length(unique(x$state))/2))
      )
    on.exit(graphics::par(oldpar))
    
    for (what in unique(x$state)) {

      graphics::boxplot(
        counts ~ date,
        data = x[x$state == what,,drop=FALSE],
        main = what,
        xlab = "Date",
        ylab = "Counts",
        border = "black",
        las = 2
        )
        
    }
    
  } else {
    
    plot.epiworld_multiple_save_reproductive_number(x, ...)

  }
  
  
}

#' @export
plot.epiworld_multiple_save_reproductive_number <- function(x, y = NULL, ...) {

  # Identifying sims
  sims <- sort(unique(x[["sim_num"]]))
  
  totals <- NULL
  for (s in sims) {
    
    # Subsetting the data
    x_tmp <- x[x[["sim_num"]] == s,, drop = FALSE]
    
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

#' @export
#' @rdname run_multiple
#' @aliases epiworld_saver
make_saver <- function(
    ...,
    fn = ""
    ) {
  
  what <- list(...)
  
  # Any missmatch?
  available <- c(
    "total_hist",
    "variant_info",
    "variant_hist",
    "tool_info",
    "tool_hist",
    "transmission",
    "transition",
    "reproductive",
    "generation"
  )
  
  not_in_available <- which(!(what %in% available))
  if (length(not_in_available)) {
    stop(
      "The following elements in -what- are not supported: \"",
      paste(what[not_in_available], collapse = "\" , \""),
      "\""
      )
  }
  
  what_bool <- as.list(available %in% what)
  names(what_bool) <- available
  
  # Checking the filename
  file_output <- TRUE
  if (fn == "") {
    fn <- file.path(tempdir(), "%05lu-episimulation.csv")
    file_output <- FALSE
  } else if (!dir.exists(dirname(fn))) {
    stop("The directory \"", dirname(fn), "\" does not exists.")
  }
  
  what_bool$fn <- fn
  
  # Generating the saver
  structure(
    list(
      ptr         = do.call(make_saver_cpp, what_bool),
      fn          = fn,
      file_output = file_output,
      what        = available[which(available %in% what)]
      ),
    class = "epiworld_saver"
  )
  
}

#' @export
print.epiworld_saver <- function(x, ...) {
  
  cat("A saver for -run_multiple-\n")
  cat("Saves the following:", paste(x$what, sep = ", "), "\n")
  cat("To file            :", ifelse(x$file_output, "yes", "no"), "\n")
  if (x$file_output)
    cat("Saver pattern      :", x$fn)
  
  invisible(x)
}

