# devtools::install_github("UofUEpi/epiworldR")
library(epiworldR)
library(data.table)

N     <- 500
n     <- 1000
ndays <- 100

set.seed(1231)

theta <- data.table(
  preval = rbeta(N, 1, 19),        # Mean 1/(1 + 19) = 0.05
  repnum = rgamma(N, 4, 4/1.5),    # Mean 4/(4 / 1.5) = 1.5
  ptran  = rbeta(N, 19, 1),        # Mean 19/(1 + 19) = 0.95
  prec   = rbeta(N, 10, 10*2 - 10) # Mean 10 / (10 * 2 - 10) = .5
)

ans <- vector("list", N)
for (i in 1:N) {
  
  m <- theta[i,
    ModelSIRCONN(
      "mycon",
      prevalence = preval,
      reproductive_number = repnum,
      prob_transmission = ptran,
      prob_recovery = prec,
      n = n
      )
    ]
  
  
  run(m, ndays = ndays)
  ans[[i]] <- get_hist_total(m)
  
}

# Setting up the data for tensorflow. Need to figure out how we would configure
# this to store an array of shape 3 x 100 (three rows, S I R) and create the 
# convolution.

# 100 Matrices of 3x101 - SIR as rows 
matrices <- parallel::mclapply(ans, function(x) {
                susceptible <- x[x[,2] == "Susceptible",]$counts
                infected <- x[x[,2] == "Infected",]$counts
                recovered <- x[x[,2] == "Recovered",]$counts
                matrix_3x101 <- matrix(rbind(susceptible, infected, recovered), nrow = 3)
            }, mc.cores = parallel::detectCores())

# Flattening out each matrix to 1-D Arrays
# Initialize a list to store the 1-dimensional arrays
arrays_1d <- vector("list", length = N)
# Iterate through each matrix
for (i in 1:length(matrices)) {
  arrays_1d[[i]] <- as.vector(matrices[[i]])
}

# Convolutional Neural Network
library(keras)

# (N obs, rows, cols)
# Important note, it is better for the model to handle changes rather than
# total numbers. For the next step, we need to do it using % change, maybe...
arrays_1d <- array(dim = c(N, 3, 50))
for (i in seq_along(matrices))
  arrays_1d[i,,] <- t(diff(t(matrices[[i]])))[,1:50]

theta2 <- copy(theta)
theta2[, repnum := plogis(repnum)]

# Reshaping
train <- list(
  x = array_reshape(arrays_1d, dim = c(N, 3, 50)),
  y = array_reshape(as.matrix(theta2), dim = c(N, 4))
)

# Follow examples in: https://tensorflow.rstudio.com/tutorials/keras/classification

# Build the model
model <- keras_model_sequential()
model %>%
  layer_conv_2d(
    filters     = 32,
    input_shape = c(3, 50, 1),
    activation  = "relu",
    kernel_size = c(3, 5)
    ) %>%
  layer_max_pooling_2d(
    pool_size = 2,
    padding = 'same'
    ) %>%
  layer_flatten(
    input_shape = c(3, 50)
    ) %>%
  layer_dense(
    units = 4,
    activation = 'sigmoid'
    )

# Compile the model
model %>% compile(
  optimizer = 'adam',
  loss = 'mse'
)

# Running the model
model %>% fit(train$x, train$y, epochs = 50, verbose = 2)

pred <- predict(model, x = train$x)
abs(pred - as.matrix(theta2)) 