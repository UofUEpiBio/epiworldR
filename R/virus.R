#' Virus
#' @param name of the virus
#' @param post_immunity,prob_infecting,prob_recovery,prob_death to be documented
#' @examples 
#' mseirconn <- ModelSEIRCONN(
#'   name                = "COVID-19",
#'   prevalence          = 0.01, 
#'   n                   = 10000,
#'   reproductive_number = 4, 
#'   incubation_days     = 7, 
#'   prob_transmission   = 0.5,
#'   prob_recovery       = 0.99
#' )
#' 
#' delta <- virus("Delta Variant", 0, .5, .2, .01)
#' virus_set_status(delta, 1, 3, 3)
#' 
#' add_virus(mseirconn, delta, .3)
#' init(mseirconn, 100, 992)
#' mseirconn
#' 
#' run(mseirconn)
#' mseirconn
#' @export
#' 
virus <- function(
    name,
    post_immunity,
    prob_infecting,
    prob_recovery,
    prob_death
    ) {
  
  structure(
    virus_cpp(
      name,
      post_immunity,
      prob_infecting,
      prob_recovery,
      prob_death
      ),
    class = "epiworld_virus"
  )
  
}

#' @export
#' @rdname virus
add_virus <- function(m, v, prevalence) {
  
  if (inherits(m, "epiworld_seirconn"))
    add_virus_cpp(m, v, prevalence)
  else
    stop("No method for object of class ", class(m))
  
  invisible(m)
}

#' @export
#' @rdname virus
add_virus_n <- function(m, v, prevalence) {
  
  if (inherits(m, "epiworld_seirconn"))
    add_virus_n_cpp(m, v, prevalence)
  else
    stop("No method for object of class ", class(m))
  
  invisible(m)
}

#' @export
#' @rdname virus
virus_set_status <- function(v, init, end, removed) {
  
  virus_set_status_cpp(v, init, end, removed)
  invisible(v)
  
}
