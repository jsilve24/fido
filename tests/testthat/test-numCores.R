set.seed(1)
library(fido)
sim <- pibble_sim()
fit <- optimPibbleCollapsed(sim$Y, sim$upsilon, sim$Theta%*%sim$X, sim$KInv, 
                            sim$AInv, random_pibble_init(sim$Y))  


test_that("Number of Cores Test 1", {
  fit2 <- uncollapsePibble(fit$Samples, sim$X, sim$Theta, 
                           sim$Gamma, sim$Xi, sim$upsilon, 
                           seed=2849, ncores = 1)
  expect_equal(fit2$NoCores, 1)
})
