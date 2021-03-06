context("printAovMeans")

test_that("printAovMeans", {

  # create dataframe
  dat <- createDF(nVP = 6, nTrl = 1,
                  design = list("Comp" = c("comp", "incomp")))
  dat <- addDataDF(dat, RT = list("Comp comp"   = c(500, 150, 100),
                                  "Comp incomp" = c(520, 150, 100)))

  # base R aov
  aovRT <- aov(RT ~ Comp + Error(VP/Comp), dat)
  aovRT <- aovTable(aovRT)

  testthat::expect_error(printAovMeans(aovRT), NA)

  # ezANOVA
  aovRT <- ez::ezANOVA(dat, dv = .(RT), wid = .(VP), within = .(Comp),
                       return_aov = TRUE, detailed = TRUE)
  aovRT <- aovTable(aovRT)
  testthat::expect_error(printAovMeans(aovRT), NA)

  aovRT <- ez::ezANOVA(dat, dv = .(RT), wid = .(VP), within = .(Comp),
                       return_aov = FALSE, detailed = TRUE)
  testthat::expect_error(printAovMeans(aovRT))

  aovRT <- ez::ezANOVA(dat, dv = .(RT), wid = .(VP), within = .(Comp),
                       return_aov = TRUE, detailed = TRUE)
  aovRT <- aovTable(aovRT)
  testthat::expect_error(printAovMeans(aovRT, aovRT, digits = c(2, 2)), NA)
  testthat::expect_error(printAovMeans(aovRT, digits = c(2, 2)))
  testthat::expect_error(printAovMeans(aovRT, dv = c("ms1", "ms2")))

})
