test_that("summary works", {
  x <- rnorm(100)
  expect_error(expect_error(summary(x))) # expect no error!
})

test_that("print works", {
  x <- rnorm(100)
  expect_error(expect_error(print(x))) # expect no error!
})

test_that("plot works", {
  x <- rnorm(100)
  y <- rnorm(100)
  expect_error(expect_error(plot(x,y))) # expect no error!
  
  set.seed(1)
  sim <- pibble_sim(N=10, D=4, Q=3)
  fit <- pibble(sim$Y, sim$X)
  expect_error(expect_error(plot(fit, par="Lambda"))) # expect no error!
})
