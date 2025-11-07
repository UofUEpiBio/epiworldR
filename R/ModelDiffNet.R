#' Network Diffusion Model
#'
#' The network diffusion model is a simple model that assumes that
#' the probability of adoption of a behavior is proportional to the
#' number of adopters in the network.
#'
#' @export
#' @concept general-models
#' @param name Name of the model.
#' @param prevalence Prevalence of the disease.
#' @param prob_adopt Probability of adoption.
#' @param normalize_exposure Normalize exposure.
#' @param data Data.
#' @param data_cols Data columns.
#' @param params Parameters.
#' @return An object of class [epiworld_diffnet] and [epiworld_model].
#' @family Models
#' @details
#' Different from common epidemiological models, the network diffusion model
#' assumes that the probability of adoption of a behavior is proportional to the
#' number of adopters in the network. The model is defined by the following
#' equations:
#' \deqn{
#' P(adopt) = \mbox{Logit}^{-1}(prob\_adopt + params * data + exposure)
#' }
#' Where exposure is the number of adopters in the agent's network.
#'
#' Another important difference is that the transmission network is not
#' necesary useful since adoption in this model is not from a particular
#' neighbor.
#'
#' @examples
#' set.seed(2223)
#' n <- 10000
#'
#' # Generating synthetic data on a matrix with 2 columns.
#' X <- cbind(
#'   age = sample(1:100, n, replace = TRUE),
#'   female = sample.int(2, n, replace = TRUE) - 1
#' )
#'
#' adopt_chatgpt <- ModelDiffNet(
#'   "ChatGPT",
#'   prevalence = .01,
#'   prob_adopt = .1,
#'   data       = X,
#'   params     = c(1, 4)
#' )
#'
#' # Simulating a population from smallworld
#' agents_smallworld(adopt_chatgpt, n, 8, FALSE, .01)
#'
#' # Running the model for 50 steps
#' run(adopt_chatgpt, 50)
#'
#' # Plotting the model
#' plot(adopt_chatgpt)
#' @aliases epiworld_diffnet
ModelDiffNet <- function(
  name,
  prevalence,
  prob_adopt,
  normalize_exposure = TRUE,
  data               = matrix(nrow = 0, ncol = 0),
  data_cols          = 1L:ncol(data),
  params             = vector("double")
) {
  # Check input params
  stopifnot_string(name)
  stopifnot_double(prevalence)
  stopifnot_double(prob_adopt)
  stopifnot_bool(normalize_exposure)
  stopifany_na(data)
  stopifnot_numvector(data_cols)
  stopifnot_numvector(params)

  if (length(data) == 0L)
    data_cols <- vector("integer")
  else
    data_cols <- as.integer(data_cols) - 1L

  structure(
    ModelDiffNet_cpp(
      name,
      prevalence,
      prob_adopt,
      normalize_exposure,
      data,
      ncol(data),
      data_cols,
      params
    ),
    class = c("epiworld_diffnet", "epiworld_model")
  )

}
