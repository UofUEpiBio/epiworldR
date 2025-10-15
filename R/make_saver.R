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
#' - `virus_info` Information about `viruses`.
#' - `virus_hist` Changes in `viruses`.
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
#'   transmission_rate = 0.9, recovery_rate = 0.1
#' )
#'
#' # Generating a saver
#' saver <- make_saver("total_hist", "reproductive")
#'
#' # Running and printing
#' run_multiple(model_sir, ndays = 100, nsims = 50, saver = saver, nthreads = 2)
#'
#' # Retrieving the results
#' ans <- run_multiple_get_results(model_sir, nthreads = 2)
#'
#' head(ans$total_hist)
#' head(ans$reproductive)
#'
#' # Plotting
#' multi_sir <- ans$total_hist
#' multi_sir <- multi_sir[multi_sir$date <= 20, ]
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

  # If the saver is greater than 1, then
  # we need to delete the files from the previous run
  fnames <- list.files(
    path = dirname(saver$fn),
    full.names = TRUE
  )

  if (length(fnames)) {
    unlink(fnames, expand = FALSE)
  }

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
#' @param nthreads Integer. Number of threads (passed to [parallel::makeCluster()]).
#' @param freader A function to read the files. If `NULL` (default,) uses
#' `utils::read.table`.
#' @param ... Additional arguments passed to `freader`.
#' @details
#' An alternative to using the default `utils::read.table` function is to use
#' `data.table::fread` from the `data.table` package. This can be done by
#' specifying `freader = data.table::fread` and passing additional arguments
#' (e.g., `nThread = 2L`) via `...`. This can significantly speed up the
#' reading process, especially for large datasets.
#'
#' If the model does not, for example, tools, then the corresponding data frame
#' will be empty (0 rows). A warning will be issued in this case when trying
#' to retrieve or plot the results.
#' @returns
#' - The `run_multiple_get_results` function returns a named list with the
#' data specified by `make_saver`.
#' @importFrom utils read.table
#' @importFrom parallel parLapply makeCluster stopCluster detectCores
run_multiple_get_results <- function(
    m,
    nthreads = min(2L, parallel::detectCores()),
    freader = NULL,
    ...
    ) {

  if (!inherits(m, "epiworld_model"))
    stop("-m- must be of class `epiworld_model`.")

  # Get the filepath
  saver <- attr(m, "saver")

  if (!length(saver))
    stop("No -saver- found. -run_multiple_get_results- can only be used after using -run_multiple-.")

  output <- vector("list", length(saver$what))
  names(output) <- saver$what

  cl <- if (nthreads > 1L)
    parallel::makeCluster(nthreads)
  else NULL

  on.exit(if (length(cl)) parallel::stopCluster(cl))

  for (i in saver$what) {
    # Listing the files
    fnames <- list.files(
      path       = dirname(saver$fn),
      pattern    = sprintf("%s\\.csv", i),
      full.names = TRUE
    )

    reader_fun <- if (length(freader))
      function(x, ...) freader(x, ...)
    else
      function(x, ...) utils::read.table(
        x,
        sep = " ", header = TRUE, comment.char = ""
      )

    # Reading the files
    output[[i]] <- if (length(cl))
      parallel::parLapply(cl, fnames, reader_fun, ...)
    else
      lapply(fnames, reader_fun, ...)

    # Getting number of simulation
    output[[i]] <- lapply(seq_along(fnames), \(j) {
      # It doesn't matter if the file is empty
      # but we can't cbind if it is empty
      if (nrow(output[[i]][[j]]) == 0) {
        cnames <- c("sim_num", colnames(output[[i]][[j]]))

        return(
          data.frame(
            matrix(
              NA,
              nrow = 0,
              ncol = length(cnames),
              dimnames = list(NULL, cnames)
            )
          )
        )

      }

      cbind(sim_num = j, output[[i]][[j]])

    })

    # Putting all together
    output[[i]] <- do.call(rbind, output[[i]])

    # If there are no observations, then
    err_msg <- tryCatch(
      {
        class(output[[i]]) <- c("epiworld_multiple_save_i", class(output[[i]]))
      },
      error = function(e) e
    )

    if (inherits(err_msg, "error")) {

      warning(
        "When retrieving the saved results, for the case of ",
        i, ", there were no observations."
      )

      class(output[[i]]) <- c(
        "epiworld_multiple_save_i",
        class(output[[i]])
      )

    }

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

  if (nrow(x) == 0) {
    warning(
      "When plotting the saved results, for the case of ",
      what, ", there were no observations."
    )
    return(NULL)
  }

  # If it is not reproductive number, then...
  if (what != "reproductive") {

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
    "virus_info",
    "virus_hist",
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

  # Using tempfile to generate directories
  id <- basename(tempfile("epiworldR-"))

  if (fn == "") {

    fp <- file.path(tempdir(), id)
    dir.create(fp)
    fn <- file.path(fp, "%05lu-episim")
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
      what        = available[which(available %in% what)],
      id          = id
    ),
    class = "epiworld_saver"
  )

}

#' @export
print.epiworld_saver <- function(x, ...) {

  cat("A saver for -run_multiple-\n")
  cat("Saves the following :", paste(x$what, sep = ", "), "\n")
  cat("To file             :", ifelse(x$file_output, "yes", "no"), "\n")
  cat("Saver folder        :", dirname(x$fn), "\n")
  cat("Saver id            :", x$id, "\n")

  invisible(x)

}
