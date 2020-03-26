N <- 10 #number of obs
P <- 100 #number of variables
data <- data.frame(matrix(rnorm(N*(P+1)), nrow = N, ncol = P+1))

set.seed(100)
sln <- FSA(formula = "X101~1", data = data, cores = 1, m = 2,
           interactions = FALSE, criterion = AIC, minmax = "min",
           numrs = 10)
sln
# formula criterion times
# Original Fit     X101 ~ 1 22.670894    NA
# FS1          X101~X32+X70  7.969452     3
# FS2          X101~X40+X84  7.747500     2
# FS3          X101~X55+X73 10.156049     5

set.seed(100)
sln2 <- FSA2(formula = "X101~1", data = data, cores = 1, m = 2,
           interactions = FALSE, criterion = AIC, minmax = "min",
           numrs = 10)
sln2
> sln2
# formula criterion times
# Original Fit     X101 ~ 1 22.670894    NA
# FS1          X101~X32+X70  7.969452     3
# FS2          X101~X40+X84  7.747500     2
# FS3          X101~X55+X73 10.156049     5