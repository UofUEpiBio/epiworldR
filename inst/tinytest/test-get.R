
model <- ModelSIR(
    name                = "COVID-19",
    n                   = 10000,
    prevalence          = 0.01,
    transmission_prob   = 0.4,
    recovery_rate       = 0.95
)
add_param(model, "Contact rate", 5)
add_layer(
  model,
  layer_mixing("everyone", groups = rep(0, 10000), contact_matrix = 1) |>
    set_contact_rate("Contact rate")
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