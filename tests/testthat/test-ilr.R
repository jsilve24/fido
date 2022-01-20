test_that("ilr plus summary works", {
  library(MicrobeDS)
  library(phyloseq)
  library(dplyr)
  
  Y = pibble_sim()
  priors <- pibble(NULL, Y$X, Y$upsilon, Y$Theta, Y$Gamma, Y$Xi)
  priors <- to_ilr(priors)
  expect_error(expect_error(summary(priors, pars="Lambda"))) # expect no error!
})

test_that("clr plus summary works", {
  library(MicrobeDS)
  library(phyloseq)
  library(dplyr)
  
  Y = pibble_sim()
  priors <- pibble(NULL, Y$X, Y$upsilon, Y$Theta, Y$Gamma, Y$Xi)
  priors <- to_clr(priors)
  expect_error(expect_error(summary(priors, pars="Lambda"))) # expect no error!
})