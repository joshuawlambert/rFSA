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
  formula = "X101~1", data = data, cores = 1, m = 2,
  interactions = F, criterion = AIC, minmax = "min",
  numrs = numrs)
summary(res)
test_that("multiplication works", {
  expect_equivalent(as.formula(res$table$formula), as.formula("X101~X7+X83"))
  expect_equal(res$table$criterion,266.028,
               tolerance=1e-2,check.attributes=F)
})
