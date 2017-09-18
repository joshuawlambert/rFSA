#' FSA: Feasible Solution Algorithm
#'
#' @description A function using a Feasible Solution Algorithm to find a set of feasible solutions for a statistical model of a specific form that could include mth-order interactions (Note that these solutions are optimal in the sense that no one swap to any of the variables will increase the criterion function.)
#' 
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. 
#' @param data a data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
#' @param fitfunc the method that should be used to fit the model. For Example: lm, glm, or other methods that rely on formula, data, and other inputs.
#' @param fixvar variable(s) to fix in the model. Usually a covariate that should always be included (Example: Age, Sex). Will still consider it with interactions. Default is NULL.
#' @param quad Include quadratic terms or not. Logical.
#' @param m order of terms to include. If interactions is set to TRUE then m is the order of interactions to be considered. For Subset selection (interaction=F), m is the size of the subset to examine. Defaults to 2.  
#' @param numrs number of random starts to perform.
#' @param cores number of cores to use while running. Note: Windows can only use 1 core. See mclapply for details. If function detects a Windows user it will automatically set cores=1.
#' @param interactions whether to include interactions in model. Defaults to TRUE.
#' @param criterion which criterion function to either maximize or minimize. For linear models one can use: r.squared, adj.r.squared, cv5.lmFSA (5 Fold Cross Validation error), cv10.lmFSA (10 Fold Cross Validation error), apress (Allen's Press Statistic), int.p.val (Interaction P-value), AIC, BIC.
#' @param minmax whether to minimize or maximize the criterion function
#' @param checkfeas vector of variables that could be a feasible solution. These variables will be used as the last random start.
#' @param var4int specification of which variables to check for marginal feasiblilty. Default is NULL
#' @param min.nonmissing the combination of predictors will be ignored unless this many of observations are not missing
#' @param ... other arguments passed to fitfunc.
#'
#' @import hashmap
#' @importFrom parallel mclapply
#' @return matrix of results
#' @export
#'
#' @examples
#' 
#' N <- 10 #number of obs
#' P <- 100 #number of variables
#' data <- data.frame(matrix(rnorm(N*(P+1)), nrow = N, ncol = P+1))
#' 
#' sln <- FSA(formula = "X101~1", data = data, cores = 1, m = 2,
#' interactions = F, criterion = AIC, minmax = "min",
#' numrs = 10)
#' sln
#'
#' @describeIn FSA find best set of variables for statistical models
FSA <- function(formula, data, fitfunc=lm, fixvar = NULL, quad = FALSE,
                m = 2, numrs = 1, cores=1, interactions = T,
                criterion = AIC, minmax="min", checkfeas=NULL, var4int=NULL,
                min.nonmissing=1,...)
{
  formula <- as.formula(formula)
  data <- data.frame(data)

  if(.Platform$OS.type != "unix") cores = 1

  if(!(is.atomic(min.nonmissing) & length(min.nonmissing)==1)) {
    stop("min.nonmissing should be a scalar.")
  }
  
  yname <- all.vars(formula)[1]
  allname <- colnames(data)
  stopifnot(all(c(yname,fixvar) %in% allname))
  P <- length(allname)-1
  ypos <- which(allname == yname)
  xpos <- setdiff(1:(P+1), ypos)
  criterion <- criterion
  which.best <- switch(tolower(minmax), min=which.min, max=which.max)
  bad.cri <- switch(tolower(minmax), min=Inf, max=-Inf)

  ##Generate random starting positions
  #if checkfeas != NULL and length(checkfeas)==m then put the the check feas in the last position of starts
  if (is.null(checkfeas)) {
    starts <- replicate(n=numrs, expr=pos2key(sort(sample(xpos, m, replace = F))))
  } else if (length(checkfeas)!=m) {
    stop("sorry, the number of variables in checkfeas is not equal to m. Please try again.")
  } else {
    starts <- replicate(n=numrs, expr=pos2key(sort(sample(xpos, m, replace = F))))
    starts[length(starts)]<-pos2key(which(colnames(data) %in% checkfeas))
  }
  
  ## cur.key stores the keys of current positions
  ## during optimization.
  cur.key <- starts

  
  form.str <- function(val)
  {
    str = paste0(yname,"~")
    if(!is.null(fixvar)) str = paste0(paste0(str,fixvar, collapse = "+"),"+")
    str=paste0(str,paste0(allname[val], collapse=ifelse(isTRUE(interactions),"*","+")))
  }
  
  form <- function(val){as.formula(form.str(val))}

  ## Initilize a hash table to store criterion for the computed
  ## combinations. The keys used to index criterions are
  ## produced by pos2key, and could be decoded by key2pos
  Cri <- hashmap("",1)
  Cri$erase("")

  info <- data.frame(
    start=starts, current=starts, solution=NA, iteration=0,
   check=0,stringsAsFactors = F)
  unsolved.mask <- is.na(info$solution)
  while(any(unsolved.mask))
  {
    info$iteration[unsolved.mask] <- info$iteration[unsolved.mask] + 1
    unsolved.cur <- info$current[unsolved.mask]

    steps <- lapply(
      unsolved.cur,
      FUN = function(x)
      {
        tmp <- swaps(cur = key2pos(x), n=P+1, quad = quad,yindex = ypos)
        if (!is.null(var4int)) {
          tmp <- tmp[,which(apply(tmp==which(allname==var4int), MARGIN=2, FUN=any))]
        }
        apply(tmp, MARGIN=2, FUN=pos2key)
      }
    )
    names(steps) <- unsolved.cur
    info$check[unsolved.mask] <- info$check[unsolved.mask] + sapply(steps,length)

    ## Calculate criterion for each next step position
    ## Basically, we will check the criterion hash table
    ## If a combination's criterion already exists in
    ## the hash table, we simply use it. Otherwise, we
    ## will calculate the criterion and insert it into
    ## the hash table
    steps.noCri <- unique(unlist(steps))
    steps.noCri <- steps.noCri[!Cri$has_key(steps.noCri)];
    if (length(steps.noCri) > 0 )
    {
      new.Cri <- unlist(mclapply(
        steps.noCri, mc.cores=cores,
        FUN = function(key)
        {
          pos <- key2pos(key)
          ## check if there are too many NAs
          if (sum(!apply(is.na(data[,c(pos, ypos)]), 1, any)) < min.nonmissing) {
            bad.cri
          } else {
            tryCatch(criterion(fitfunc(formula=form(pos), data = data,...)),
                     error=function(cond) {bad.cri})
          }
        }))
      Cri[[steps.noCri]] <- new.Cri
    }
    
    ##Find the best next position for each current position
    unsolved.next <- unlist(mclapply(
      unsolved.cur, mc.cores=cores,
      FUN = function(key)
      {
        step <- steps[[key]]
        step[which.best(Cri[[step]])]
      }
    ))
    stopifnot(length(unsolved.cur)==length(unsolved.next))

    ##Check if any solutions are found
    mask <- (unsolved.cur == unsolved.next |
               unsolved.next %in% info$solution)
    #mask <- unsolved.cur == unsolved.next
    
    cur.sln <- rep(NA, length(unsolved.cur))
    cur.sln[mask] <- unsolved.next[mask]

    ##Update settings and iterate
    info$solution[unsolved.mask] <- cur.sln
    info$current[unsolved.mask] <- unsolved.next
    unsolved.mask <- is.na(info$solution)
  }


  ##************************************************************
  ## format outputs
  ##************************************************************
  originalfit <- tryCatch(fitfunc(formula=formula, data=data),
                          error=function(e){NULL})
  call <- mget(names(formals()),sys.frame(sys.nframe()))

  solutions <- list()
  for (k in 1:m) {
    solutions[[paste0("start.",k)]] <- sapply(
      info$start,
      FUN=function(key){allname[key2pos(key)[k]]})
  }
  for (k in 1:m) {
    solutions[[paste0("best.",k)]] <- sapply(
      info$solution,
      FUN=function(key){allname[key2pos(key)[k]]})
  }
  solutions$criterion <- info$criterion
  solutions$swaps <- info$iteration
  solutions <- data.frame(solutions, stringsAsFactors=F)
  rownames(solutions) <- NULL
  
  sln.summary <- table(info$solution)
  sln.keys <- names(sln.summary)
  table <- list()
  table$formula <- sapply(sln.keys, FUN=function(key){form.str(key2pos(key))})
  for (k in 1:m) {
    table[[paste0("Var",k)]] <- sapply(
      sln.keys,
      FUN=function(key){allname[key2pos(key)[k]]})
  }
  table$criterion <- Cri[[sln.keys]]
  table$times <- as.numeric(sln.summary)
  table <- data.frame(table, stringsAsFactors = F)
  rownames(table) <- NULL

  efficiency <- paste(
    "You did",sum(Cri$size()), "model fittings and",
    sum(info$check), "model checks compared to",
    choose(n=P,k=m),"fittings and",choose(n=P,k=m),
    "checks you would have done with exhaustive search.")

  res <- list(originalfit=originalfit, call=call,
              solutions=solutions, table=table,
              efficiency=efficiency, info=info,
              nfits=Cri$size())
  class(res) <- "FSA"

  return(res)
}

#' @export
#' @describeIn FSA alias for \code{FSA(fitfunc=lm,...)}
lmFSA <- function(...) {FSA(fitfunc = lm, ...)}

#' @export
#' @describeIn FSA alias for \code{FSA(fitfunc=glm,...)}
glmFSA <- function(...) {FSA(fitfunc = glm, ...)}
