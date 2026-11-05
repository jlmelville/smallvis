# Rounded reference costs and finite differences use absolute numerical error.
# scale = 1 preserves the comparison used by these tests before edition 3.
expect_equal_abs <- function(object, expected, tolerance, info = NULL) {
  comparison <- all.equal(object, expected, tolerance = tolerance, scale = 1)
  testthat::expect_true(
    isTRUE(comparison),
    info = paste(c(info, comparison), collapse = "\n")
  )
}
