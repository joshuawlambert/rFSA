context("FSA")

## TODO: Rename context
## TODO: Add more tests

set.seed(1000)
numrs <- 10 #number of start point
N <- 100 #number of obs
P <- 100 #number of variables
data <- data.frame(matrix(rnorm(N*(P+1)), nrow = N, ncol = P+1))


for (cores in c(1, parallel::detectCores())) {
  test_that(paste0("random data, fixvar=NULL, cores=", cores), {
    skip_on_cran()

    set.seed(1000)
    res <- FSA(
    formula = "X101~1", data = data, cores = cores, m = 2,
    interactions = F, criterion = AIC, minmax = "min",
    numrs = numrs, return.models=FALSE)

    expect_equivalent(as.formula(res$table$formula), as.formula("X101~X7+X83"))
    expect_equal(res$table$criterion,266.028,
                 tolerance=1e-2,check.attributes=F)

    set.seed(1000)
    res <- FSA(
      formula = "X101~1", data = data, cores = cores, m = 2,
      interactions = F, criterion = AIC, minmax = "min",
      numrs = numrs, return.models=TRUE)

    expect_equivalent(as.formula(res$table$formula), as.formula("X101~X7+X83"))
    expect_equal(res$table$criterion,266.028,
                 tolerance=1e-2,check.attributes=F)
    expect_true("swapped.to.model" %in% names(res$solutions))
    expect_true("checked.model" %in% names(res$solutions))
    
  })
}

for (cores in c(1, parallel::detectCores())) {
  test_that(paste0("mtcars data, fixvar=hp, cores=", cores),{
    skip_on_cran()
    
    set.seed(1000)
    res <- FSA(fitfunc=lm,formula="mpg~hp+wt", data=mtcars, fixvar="hp",
               quad=F, m=2, numrs=10, cores=1, criterion = r.squared, minmax="max",
               return.models=FALSE)
    expect_equivalent(as.formula(res$table$formula), as.formula("mpg~hp+hp*wt"))
    expect_equal(res$table$criterion, 0.8847637, tolerance=1e-3,check.attributes=F)

    set.seed(1000)
    res <- FSA(fitfunc=lm,formula="mpg~hp+wt", data=mtcars, fixvar="hp",
               quad=F, m=2, numrs=10, cores=1, criterion = r.squared, minmax="max",
               return.models=TRUE)
    expect_equivalent(as.formula(res$table$formula), as.formula("mpg~hp+hp*wt"))
    expect_equal(res$table$criterion, 0.8847637, tolerance=1e-3,check.attributes=F)
    expect_true("swapped.to.model" %in% names(res$solutions))
    expect_true("checked.model" %in% names(res$solutions))
  })
}
