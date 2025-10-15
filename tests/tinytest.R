
if (requireNamespace("tinytest", quietly = TRUE)) {

  test_tmat <- function(tmat) {
    # Check for out of bounds values
    expect_false(any(tmat < 0))
    expect_false(any(tmat > 1))
  }

  tinytest::test_package("epiworldR")
}
