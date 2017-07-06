#examples

###************************************************************
### performance benchmark
###************************************************************
numrs <- 100 #number of start point
N <- 100 #number of obs
P <- 100 #number of variables

data <- data.frame(matrix(rnorm(N*(P+1)), nrow = N, ncol = P+1))


t1 <- Sys.time()
sln <- FSA(formula = "X101~1", data = data, cores = 1, m = 2,
           interactions = F, criterion = AIC, minmax = "min",
           numrs = numrs,usehist=T)
Sys.time()-t1

t1 <- Sys.time()
sln <- FSA(formula = "X101~1", data = data, cores = 1, m = 2,
           interactions = F, criterion = AIC, minmax = "min",
           numrs = numrs,usehist=F)
Sys.time()-t1


t1 <- Sys.time()
fit <- lmFSA(formula = "X101~1", data = data, cores = 1, m = 2,
             interactions = F, criterion = AIC, minmax = "min",
             numrs = numrs)
Sys.time()-t1
