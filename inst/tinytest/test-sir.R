# Adding a Small world population without queuing ------------------------------
sir_0 <- ModelSIR(
  name = "COVID-19",
  prevalence = .01,
  infectiousness = .9,
  recovery = .3
)

agents_smallworld(
  sir_0,
  n = 50000,
  k = 5,
  d = FALSE,
  p = .01
)

# Initializing 
queuing_off(sir_0)
init(sir_0, days = 50, seed = 1912)

# Running and printing
run(sir_0)

tmat_0 <- get_transition_probability(sir_0)

# Creating a SIR model with queuing --------------------------------------------
sir_1 <- ModelSIR(
  name = "COVID-19",
  prevalence = .01,
  infectiousness = .9,
  recovery = .3
)

# Adding a Small world population 
agents_smallworld(
  sir_1,
  n = 50000,
  k = 5,
  d = FALSE,
  p = .01
)

# Initializing 
init(sir_1, days = 50, seed = 1912)

# Running and printing
run(sir_1)

tmat_1 <- get_transition_probability(sir_1)

# Expected
tmat_expected <- structure(c(0.90395575761795, 0, 0, 0.0959785506129265, 0.707067608833313, 
                             0, 6.55534968245775e-05, 0.292932420969009, 1), dim = c(3L, 3L
                             ), dimnames = list(c("Susceptible", "Infected", "Recovered"), 
                                                c("Susceptible", "Infected", "Recovered")))

expect_equivalent(tmat_0, tmat_1)
expect_equivalent(tmat_0, tmat_expected)
expect_equivalent(tmat_1, tmat_expected)


