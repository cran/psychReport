context("mathString")

test_that("mathString", {

  result <- mathString("2", "1", operation = "add", unit = "ms")
  expect_equal(result, "3 ms")

  result <- mathString("2 bananas", "1 apple", operation = "subtract", unit = "ms")
  expect_equal(result, "1 ms")

  result <- mathString("2 bananas", "1 apple", operation = "subtract", unit = "mv")
  expect_equal(result, "1 $\\\\mu$V")

  result <- mathString("922.2567", "621.2134", operation = "add", numDigits = 0, unit = "ms")
  expect_equal(result, "1543 ms")

  result <- mathString("9.27", "6.24", operation = "subtract", numDigits = 2, unit = "%")
  expect_equal(result, "3.03 \\\\%")

})
