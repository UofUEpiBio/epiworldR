#' Run multiple simulations at once
#' 
#' @param m,ndays,seed See [run].
#' @param saver An object of class [epiworld_saver].
#' @param nsims Integer. Number of replicats
#' @param what Character vector with the things to save (see details)
#' @param fn A file name pattern.
#' @param nthreads Currently ignored
#' 
#' @details
#' Currently, the following elements can be saved:
#' 
#' - `total_hist` 
#' - `variant_info` 
#' - `variant_hist` 
#' - `tool_info` 
#' - `tool_hist` 
#' - `transmission` 
#' - `transition` 
#' - `reproductive` 
#' - `generation` 
#' 
#' @returns
#' In the case of `make_saver`, an list of class `epiworld_saver`.
#' The function `run_multiple` returns a list of the same class as `m`; usually
#' an `epiworld_modwel` object..
#' @export
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
    seed = sample.int(1e4, 1),
    saver = make_saver(),
    reset = TRUE,
    verbose = TRUE,
    nthreads = 1L
) {
  
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
  
  invisible(m)
  
}

#' @export
#' @rdname run_multiple
read_multiple_runs <- function(saver) {
  
  if (!inherits(saver, "epiworld_saver"))
    stop("-saver- must be of class `epiworld_saver`.")
  
  # Get the filepath
  pattern <- gsub("%[0-9]*lu", "*", saver$fn)
  
  output <- vector("list", length(saver$what))
  names(output) <- saver$what
  
  for (i in saver$what) {
    
    # Listing the files
    fnames <- list.files(
      path    = dirname(pattern),
      pattern = sprintf(basename(pattern)),
      full.names = TRUE
    )
    
    # Reading the files
    output[[i]] <- lapply(fnames, read.table, sep = " ", header = TRUE)
    
    # Getting number of simulation
    output[[i]] <- lapply(seq_along(fnames), function(j) {
      cbind(sim_num = j, output[[i]][[j]])
    })
    
    # Putting all together
    output[[i]] <- do.call(rbind, output[[i]])
    
  }
  
  output
  
}

#' @export
#' @rdname run_multiple
make_saver <- function(
    what = c("total_hist"),
    fn = ""
    ) {
  
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

