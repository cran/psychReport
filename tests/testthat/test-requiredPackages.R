context("requiredPackages")

test_that("requiredPackages", {

  expect_error(requiredPackages(c("ez", "reshape2")), NA)
  expect_error(requiredPackages(c("ez", "reshape2", "xxx")))

})
