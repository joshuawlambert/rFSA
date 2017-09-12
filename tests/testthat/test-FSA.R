context("FSA")

## TODO: Rename context
## TODO: Add more tests

set.seed(1000)
numrs <- 10 #number of start point
N <- 100 #number of obs
P <- 100 #number of variables
data <- data.frame(matrix(rnorm(N*(P+1)), nrow = N, ncol = P+1))
set.seed(1000)
res <- FSA(
  formula = "X101~1", data = data, cores = 10, m = 2,
  interactions = F, criterion = AIC, minmax = "min",
  numrs = numrs)
test_that("multiplication works", {
  expect_equivalent(as.formula(res$table$formula), as.formula("X101~X7+X83"))
  expect_equal(res$table$criterion,266.028,
               tolerance=1e-2,check.attributes=F)
})

set.seed(1000)
res <- FSA(fitfunc=lm,formula="mpg~hp+wt", data=mtcars, fixvar="hp",
           quad=F, m=2, numrs=10, cores=1, criterion = r.squared, minmax="max")
test_that("fixvar",{
  expect_equivalent(as.formula(res$table$formula), as.formula("mpg~hp+hp*wt"))
  expect_equal(res$table$criterion, 0.8847637, tolerance=1e-3,check.attributes=F)
})
