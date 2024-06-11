
model <- ModelSIRCONN(
    name                = "COVID-19",
    n                   = 10000,
    prevalence          = 0.01,
    contact_rate        = 5,
    transmission_rate   = 0.4,
    recovery_rate       = 0.95
)

set.seed(1231)

# All the following functions should return a
# data.frame of length zero. Use tinytest::expect_equal()
# to check that.
expect_equal(nrow(get_hist_tool(model)), 0L)
expect_equal(nrow(get_hist_total(model)), 0L)
expect_equal(nrow(get_hist_virus(model)), 0L)
expect_equal(nrow(get_reproductive_number(model)), 0L)
expect_equal(dim(get_transition_probability(model)), c(3, 3))
expect_equal(nrow(get_hist_transition_matrix(model)), 0L)

rm(model)