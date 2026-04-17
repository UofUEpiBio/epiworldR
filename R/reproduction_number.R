#' Compute a reproduction number from a daily contact matrix with group-specific
#' infectiousness and susceptibility
#'
#' Compute the dominant eigenvalue of a next-generation matrix built from a
#' daily contact matrix, a baseline per-contact transmission probability, an
#' average infectious period, and optional group-specific infectiousness and
#' susceptibility multipliers.
#'
#' Let \eqn{C} be a daily contact matrix where \eqn{c_{ij}} is the average number
#' of daily contacts made by one infectious individual in group \eqn{i} with
#' individuals in group \eqn{j}. If \eqn{p} is a baseline infection probability
#' per contact, \eqn{D} is the mean infectious period in days, \eqn{b_i} is the
#' infectiousness multiplier for source group \eqn{i}, and \eqn{a_j} is the
#' susceptibility multiplier for recipient group \eqn{j}, then the expected
#' next-generation matrix is
#'
#' \deqn{K = p \cdot D \cdot \mathrm{diag}(b_1, \ldots, b_n) \cdot C \cdot \mathrm{diag}(a_1, \ldots, a_n)}
#'
#' or element-wise
#'
#' \deqn{k_{ij} = p \cdot D \cdot b_i \cdot c_{ij} \cdot a_j}
#'
#' The reproduction number is the spectral radius of \eqn{K}, that is, its
#' dominant eigenvalue:
#'
#' \deqn{R = \rho(K)}
#'
#' When `infectiousness` and `susceptibility` are both vectors of ones, the
#' function reduces to the homogeneous case
#'
#' \deqn{R = \rho(p \cdot D \cdot C)}
#'
#' Vaccination can be represented indirectly through these inputs. For example,
#' if vaccination reduces susceptibility only, a common choice is
#' \eqn{a_j = 1 - e_j v_j}, where \eqn{v_j} is coverage and \eqn{e_j} is vaccine
#' efficacy against susceptibility in group \eqn{j}.
#'
#' Because the contact matrix is assumed to be daily, `infectious_period_days`
#' should represent the mean number of days an infected individual remains
#' infectious enough to transmit. In simple memoryless models this is often
#' approximated by \eqn{D = 1 / \gamma}, where \eqn{\gamma} is a per-day recovery
#' rate in continuous time or, approximately, a daily recovery probability in a
#' discrete-time model.
#'
#' @param contact_matrix A square numeric matrix. Entry `[i, j]` must be the
#'   average number of daily contacts made by an infectious individual in group
#'   `i` with individuals in group `j`.
#' @param transmission_prob Numeric scalar in `[0, 1]`. Baseline infection
#'   probability per contact.
#' @param infectious_period_days Positive numeric scalar. Mean infectious period,
#'   in days. Since the contact matrix is daily, this scales daily contact
#'   opportunities to expected secondary infections over the infectious period.
#' @param infectiousness Optional numeric vector of length
#'   `nrow(contact_matrix)`. Entry `i` is the relative infectiousness multiplier
#'   for source group `i`. If `NULL`, a vector of ones is used.
#' @param susceptibility Optional numeric vector of length
#'   `nrow(contact_matrix)`. Entry `j` is the relative susceptibility multiplier
#'   for recipient group `j`. If `NULL`, a vector of ones is used.
#' @param check_reciprocity Logical. If `TRUE`, the function issues a warning
#'   when the contact matrix is not symmetric. This is only a diagnostic and
#'   does not alter the calculation. Note that reciprocity for empirical contact
#'   matrices is often assessed after adjusting for group sizes, so lack of
#'   symmetry is not automatically an error.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{R}{Numeric scalar. Dominant eigenvalue (spectral radius) of the next-generation matrix.}
#'   \item{type}{Character string, either `"R0"` for the homogeneous case or `"Reff"` when group-specific modifiers are used.}
#'   \item{next_generation_matrix}{Numeric matrix used in the calculation.}
#'   \item{infectiousness}{Numeric vector of source-group infectiousness multipliers.}
#'   \item{susceptibility}{Numeric vector of recipient-group susceptibility multipliers.}
#'   \item{eigenvalues}{Complex vector of eigenvalues of the next-generation matrix.}
#' }
#'
#' @details
#' This function implements a next-generation-matrix calculation. It is most
#' appropriate near the start of an outbreak, when the susceptible composition
#' and contact structure are approximately stable over the time window of
#' interest.
#'
#' The calculation assumes:
#' \itemize{
#'   \item the contact matrix is already on a per-day basis,
#'   \item `transmission_prob` is a baseline per-contact infection probability,
#'   \item `infectiousness` modifies the source side of transmission,
#'   \item `susceptibility` modifies the recipient side of transmission,
#'   \item a mean infectious period is an adequate summary of infectious duration.
#' }
#'
#' Interpretation note: if `infectiousness` and `susceptibility` are both equal
#' to one for all groups and the population is fully susceptible, the returned
#' value corresponds to \eqn{R_0}. Otherwise, it is better interpreted as an
#' initial or effective reproduction number under the supplied group-specific
#' transmission modifiers.
#'
#' For the per-contact interpretation to remain meaningful, the implied
#' group-pair-specific transmission probabilities
#'
#' \deqn{p_{ij} = p \cdot b_i \cdot a_j}
#'
#' should ideally remain in \eqn{[0, 1]}.
#'
#' @references
#' Diekmann O, Heesterbeek JAP, Roberts MG (2010).
#' "The construction of next-generation matrices for compartmental epidemic models."
#' Journal of the Royal Society Interface, 7(47), 873-885.
#' doi:10.1098/rsif.2009.0386
#'
#' Wallinga J, Teunis P, Kretzschmar M (2006).
#' "Using data on social contacts to estimate age-specific transmission parameters
#' for respiratory-spread infectious agents."
#' American Journal of Epidemiology, 164(10), 936-944.
#' doi:10.1093/aje/kwj317
#'
#' van den Driessche P (2017).
#' "Reproduction numbers of infectious disease models."
#' Infectious Disease Modelling, 2(3), 288-303.
#' doi:10.1016/j.idm.2017.06.002
#'
#' @examples
#' C <- matrix(c(
#'   8, 2, 1,
#'   3, 7, 2,
#'   1, 2, 5
#' ), nrow = 3, byrow = TRUE)
#'
#' # Homogeneous case
#' compute_reproduction_number(
#'   contact_matrix = C,
#'   transmission_prob = 0.04,
#'   infectious_period_days = 5
#' )
#'
#' # Group-specific infectiousness and susceptibility
#' compute_reproduction_number(
#'   contact_matrix = C,
#'   transmission_prob = 0.04,
#'   infectious_period_days = 5,
#'   infectiousness = c(1.0, 0.8, 1.2),
#'   susceptibility = c(0.5, 0.9, 0.7)
#' )
#'
#' @export
compute_reproduction_number <- function(
  contact_matrix,
  transmission_prob,
  infectious_period_days = 1,
  infectiousness = NULL,
  susceptibility = NULL,
  check_reciprocity = FALSE
) {

  if (!is.matrix(contact_matrix) || !is.numeric(contact_matrix)) {
    stop("'contact_matrix' must be a numeric matrix.")
  }

  if (nrow(contact_matrix) != ncol(contact_matrix)) {
    stop("'contact_matrix' must be square.")
  }

  if (any(!is.finite(contact_matrix))) {
    stop("'contact_matrix' must contain only finite values.")
  }

  if (any(contact_matrix < 0)) {
    stop("'contact_matrix' cannot contain negative values.")
  }

  n_groups <- nrow(contact_matrix)

  if (!is.numeric(transmission_prob) || length(transmission_prob) != 1L ||
    !is.finite(transmission_prob) || transmission_prob < 0 || transmission_prob > 1) {
    stop("'transmission_prob' must be a single number in [0, 1].")
  }

  if (
    !is.numeric(infectious_period_days) ||
      length(infectious_period_days) != 1L ||
      !is.finite(infectious_period_days) ||
      (infectious_period_days <= 0)
  ) {
    stop("'infectious_period_days' must be a single positive number.")
  }

  if (
    !is.logical(check_reciprocity) ||
      (length(check_reciprocity) != 1L) ||
      is.na(check_reciprocity)
  ) {
    stop("'check_reciprocity' must be TRUE or FALSE.")
  }

  if (
    check_reciprocity &&
      !isTRUE(all.equal(contact_matrix, t(contact_matrix)))
  ) {
    warning(
      paste(
        "'contact_matrix' is not symmetric.",
        "This may be fine if directionality is intentional,",
        "but empirical contact matrices are often balanced separately by group sizes."
      ),
      call. = FALSE
    )
  }

  if (is.null(infectiousness)) {
    infectiousness <- rep(1, n_groups)
  } else {
    if (
      !is.numeric(infectiousness) || length(infectiousness) != n_groups ||
        any(!is.finite(infectiousness)) || any(infectiousness < 0)
    ) {
      stop(
        "'infectiousness' must be NULL or a non-negative numeric vector ",
        "of length nrow(contact_matrix)."
      )
    }
  }

  if (is.null(susceptibility)) {
    susceptibility <- rep(1, n_groups)
  } else {
    if (
      !is.numeric(susceptibility) || length(susceptibility) != n_groups ||
        any(!is.finite(susceptibility)) || any(susceptibility < 0)
    ) {
      stop(
        "'susceptibility' must be NULL or a non-negative numeric vector ",
        "of length nrow(contact_matrix)."
      )
    }
  }

  implied_pair_prob <- transmission_prob * outer(infectiousness, susceptibility)

  if (any(implied_pair_prob > 1)) {
    warning(
      paste(
        "Some implied pair-specific per-contact transmission probabilities",
        "exceed 1. Check 'transmission_prob', 'infectiousness', and",
        "'susceptibility' for scale consistency."
      ),
      call. = FALSE
    )
  }

  next_generation_matrix <-
    transmission_prob * infectious_period_days *
      diag(infectiousness, nrow = n_groups, ncol = n_groups) %*%
        contact_matrix %*%
        diag(susceptibility, nrow = n_groups, ncol = n_groups)

  eigenvalues <- eigen(next_generation_matrix, only.values = TRUE)$values
  R_value <- max(Mod(eigenvalues))

  homogeneous_case <- all(infectiousness == 1) && all(susceptibility == 1)
  type <- if (homogeneous_case) "R0" else "Reff"

  list(
    R = as.numeric(R_value),
    type = type,
    next_generation_matrix = next_generation_matrix,
    infectiousness = infectiousness,
    susceptibility = susceptibility,
    eigenvalues = eigenvalues
  )
}
