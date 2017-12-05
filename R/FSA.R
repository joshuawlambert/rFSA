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
#' @param return.models bool value to specify whether return all the fitted models which have been checked
#' @param ... other arguments passed to fitfunc.
#'
#' @import hashmap
#' @importFrom parallel mclapply
#' @import tibble
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
#' interactions = FALSE, criterion = AIC, minmax = "min",
#' numrs = 10)
#' sln
#'
#' @describeIn FSA find best set of variables for statistical models
FSA <- function(formula, data, fitfunc=lm, fixvar = NULL, quad = FALSE,
                m = 2, numrs = 1, cores=1, interactions = T,
                criterion = AIC, minmax="min", checkfeas=NULL, var4int=NULL,
                min.nonmissing=1, return.models=FALSE, ...)
{
  call <- match.call()
  original <- list()
  original$formula <- Reduce(paste, deparse(as.formula(formula)))
  original$model <- fitfunc(formula=original$formula, data=data,...)
  original$model <- tryCatch(fitfunc(formula=original$formula, data=data, ...),
                             error=function(e){
                               warning("failed to fit the original model specified by formula")
                               NULL})


  if (!(is.function(criterion) | is.list(criterion) & all(sapply(criterion,is.function)))) {
    stop("criterion should be a function or a list of functions")
  }
  if (!(is.character(minmax) | is.list(minmax) & all(sapply(minmax,is.character)))) {
    stop("minmax should be character vector or list of character strings")
  }
  minmax <- tolower(unlist(minmax))
  if (!all(minmax %in% c("min","max"))) {
    stop("minmax should contain \"min\" or \"max\" only")
  }
  
  if (length(criterion) != length(minmax)) {
    stop("the number of criterion functions and number of minmax options does not match")
  }
  if (length(criterion) == 1) {
    criterion = list(criterion)
  }
  
  for (k in 1:length(criterion)) {
    if (is.null(original$model)) {
      original[[paste0("criterion.",k)]] <- NA
    } else {
      original[[paste0("criterion.",k)]] <- criterion[[k]](original$model)
    }
      
  }

  solutions <- NULL
  table <- NULL
  nfits <- 0
  nchecks <- 0
  for (k in 1:length(criterion)) {
    fit <- fitFSA(formula, data, fitfunc=fitfunc, fixvar=fixvar, quad=quad,
                  m=m, numrs=numrs, cores=cores, interactions=interactions,
                  criterion=criterion[[k]], minmax=minmax[k], checkfeas=checkfeas,
                  var4int=var4int, min.nonmissing=min.nonmissing,
                  return.models=return.models, ...)
    
    ## merge table
    table.0 <- fit$table
    N.table <- nrow(table.0)
    for (n in 1:N.table) {
      idx <- which(table.0$formula[n]==table$formula)
      if (length(idx)==0) {
        tmp <- table.0[n,]
        for (l in 1:length(criterion)) {
          tmp[[paste0("criterion.",l)]] <- criterion[[l]](fitfunc(tmp$formula, data=data, ...))
        }
        tmp$opt.criterion <- list(k)
        table <- rbind(table, tmp)
      } else if (length(idx)==1) {
        table$opt.criterion[[idx]] <- c(table$opt.criterion[idx],k)
        table$times[idx] <- table$times[idx] + tmp$times
      } else {
        stop("merge table: formulas are not unique")
      }
    }

    ## merge solutions
    solutions.0 <- fit$solution
    N.solutions <- nrow(solutions.0)
    for (l in 1:length(criterion)) {
      solutions.0[[paste0("criterion.",l)]] <- rep(NA, N.solutions)
      for (n in 1:N.solutions) {
        mask <- rep(TRUE, N.table)
        for (mm in 1:m) {
          mask <- mask &
            solutions.0[[paste0("best.",mm)]][n] == table.0[[paste0("Var",mm)]]
        }
        idx = which(mask)
        stopifnot(length(idx)==1)
        solutions.0[[paste0("criterion.",l)]][n]<- table[[paste0("criterion.",l)]][idx]
      }
    }
    solutions.0[["optimized.by"]] <- rep(paste0("criterion.",k), N.solutions)
    solutions <- rbind(solutions, solutions.0)
    
    nfits <- nfits + fit$nfits
    nchecks <- nchecks + fit$nchecks
  }

  P <- ncol(data)
  efficiency <- paste(
    "You did",nfits, "model fittings and",
    nchecks, "model checks compared to",
    length(criterion)*choose(n=P,k=m),"fittings and",
    length(criterion)*choose(n=P,k=m),
    "checks you would have done with exhaustive search.")

  res <- list(original=original, call=call,
              solutions=solutions, table=table,
              efficiency=efficiency, nfits=nfits, nchecks=nchecks)
  class(res) <- "FSA"
  return(res)
}


fitFSA <- function(formula, data, fitfunc=lm, fixvar = NULL, quad = FALSE,
                m = 2, numrs = 1, cores=1, interactions = TRUE,
                criterion = AIC, minmax="min", checkfeas=NULL, var4int=NULL,
                min.nonmissing=1, return.models=FALSE, ...)
{
  ##************************************************************
  ## check inputs
  ##************************************************************
  formula <- as.formula(formula)

  data <- as.data.frame(data)
  allname <- colnames(data)
  ##yname <- all.vars(formula)[1]
  if (!all(all.vars(formula) %in% c(allname, "."))) {
    stop(paste("variable",
               all.vars(formula)[!(all.vars(formula)%in%c(allname,"."))],
               "does not exist in data"))
  }
  yname <- all.vars(formula)[1]
  P <- length(allname)-1
  ypos <- which(allname == yname)
  xpos <- setdiff(1:(P+1), ypos)
  
  if (!is.function(fitfunc)) {
    stop("fitfunc should be a function")
  }

  if (!is.null(fixvar) && !is.character(fixvar)) {
    stop("fixvar should be NULL or a character vector")
  } else if (!all(fixvar %in% allname)) {
    stop(paste("fixvar", fixvar[!(fixvar %in% allname)],
               "does not exist in data"))
  }

  if (!is.logical(quad) | is.na(quad) | length(quad)!=1 ) {
    stop("quad should be TRUE or FALSE")
  }

  if (!(is.numeric(m) | is.integer(m)) | length(m)!=1) {
    stop("m should be a scalar")
  } else if (m<2) {
    stop("m should be greater than or equal 2")
  }

  if (!(is.numeric(numrs) | is.integer(numrs)) | length(numrs)!=1) {
    stop("numrs should be a scalar")
  } else if (numrs<1) {
    stop("numrs should be greater than or equal 1")
  }
    
  if (!(is.numeric(cores) | is.integer(cores)) | length(cores)!=1) {
    stop("cores should be a scalar")
  } else if (cores<1) {
    stop("cores should be greater than or equal 1")
  }
  if (.Platform$OS.type != "unix" & cores != 1) {
    warning("non-unix systems, force cores to be 1")
    cores = 1
  }


  if (!is.logical(interactions) | is.na(interactions) | length(interactions)!=1 ) {
    stop("interactions should be TRUE or FALSE")
  }

  if (!is.function(criterion)) {
    stop("criterion should be a function")
  }

  ## minmax has been checked in FSA()
  stopifnot(is.character(minmax) && length(minmax)==1 && minmax %in% c("min","max"))

  if (is.character(checkfeas)) {
    if (length(checkfeas)!=m) {
      stop("sorry, the number of variables in checkfeas is not equal to m. Please try again.")
    } else if (!all(checkfeas %in% allname)) {
      stop(paste("variable", checkfeas[!(checkfeas %in% allname)], "does not exist in data"))
    }
  } else if (!is.null(checkfeas)) {
    stop("checkfeas should be NULL or a character vector")
  }


  if (!(is.null(var4int) | is.character(var4int) & length(var4int)==1)) {
    stop("var4int should be NULL or a character scalar")
  }
  
  if (!((is.numeric(min.nonmissing) | is.integer(min.nonmissing))
    & length(min.nonmissing)==1)) {
    stop("min.nonmissing should be a scalar number.")
  }

  if (!is.logical(return.models) | is.na(return.models) | length(return.models)!=1 ) {
    stop("return.models should be TRUE or FALSE")
  }
  
  which.best <- switch(tolower(minmax), min=which.min, max=which.max)
  bad.cri <- switch(tolower(minmax), min=Inf, max=-Inf)

  ##Generate random starting positions
  #if checkfeas != NULL and length(checkfeas)==m then put the the check feas in the last position of starts
  stopifnot(is.null(checkfeas) | all(checkfeas %in% allname))
  if (is.null(checkfeas)) {
    starts <- replicate(n=numrs, expr=pos2key(sort(sample(xpos, m, replace = F))))
  } else {
    starts <- replicate(n=numrs, expr=pos2key(sort(sample(xpos, m, replace = F))))
    starts[length(starts)]<-pos2key(which(colnames(data) %in% checkfeas))
  }

  ##************************************************************
  ## optimization
  ##************************************************************

  ## cur.key stores the keys of current positions
  ## during optimization.
  cur.key <- starts

  form.str <- function(val)
  {
    str = paste0(yname,"~")
    if(!is.null(fixvar)) {
      str=paste0(str,paste0(fixvar,collapse="+"),"+")
    }
    str=paste0(str,paste0(allname[val], collapse=ifelse(isTRUE(interactions),"*","+")))
  }
  form <- function(val){as.formula(form.str(val))}

  ## Initilize a hash table to store criterion for the computed
  ## combinations. The keys used to index criterions are
  ## produced by pos2key, and could be decoded by key2pos
  Cri <- hashmap("",1)
  Cri$erase("")
  if (isTRUE(return.models)) {
    MDL <- list()
  }

  info <- tibble(
    start=starts, current=starts, solution=NA, iteration=0,
    check=0, steps=as.list(starts), history=as.list(starts))

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
    
    for (k in 1:length(unsolved.cur)) {
      idx.global <- (1:numrs)[unsolved.mask][k]
      info$steps[[idx.global]] <- c(info$steps[[idx.global]], steps[[k]])
    }

    ## Calculate criterion for each next step position
    ## Basically, we will check the criterion hash table
    ## If a combination's criterion already exists in
    ## the hash table, we simply use it. Otherwise, we
    ## will calculate the criterion and insert it into
    ## the hash table
    steps.noCri <- unique(unlist(steps))
    steps.noCri <- steps.noCri[!Cri$has_keys(steps.noCri)];

    if (length(steps.noCri) > 0 )
    {
      new.Cri <- mclapply(
        steps.noCri, mc.cores=cores,
        FUN = function(key)
        {
          pos <- key2pos(key)

          if (isTRUE(return.models)) {
            if (sum(!apply(is.na(data[,c(pos, ypos)]), 1, any)) < min.nonmissing) {
              list(criterion=bad.cri, model=NA)
            } else {
              tryCatch({
                mdl <- fitfunc(formula=form(pos), data=data,...);
                cri <- criterion(mdl);
                list(criterion=cri, model=mdl)
              },
              error=function(cond) {list(criterion=bad.cri, model=NA)})
            }
          } else {
            if (sum(!apply(is.na(data[,c(pos, ypos)]), 1, any)) < min.nonmissing) {
              bad.cri
            } else {
              tryCatch(
                criterion(fitfunc(formula=form(pos), data = data,...)),
                error=function(cond) {bad.cri})
            }
          }
        }
        )

      if (isTRUE(return.models)) {
        for (k in 1:length(steps.noCri)) {
          Cri[[steps.noCri[k]]] <- new.Cri[[k]]$criterion
          MDL[[steps.noCri[k]]] <- new.Cri[[k]]$model
        }
      } else {
        Cri[[steps.noCri]] <- unlist(new.Cri)
      }
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

    for (k in 1:length(unsolved.cur)) {
      idx.global <- (1:numrs)[unsolved.mask][k]
      info$history[[idx.global]] <- c(info$history[[idx.global]], unsolved.next[k])
    }
    
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
  

  originalfit <- tryCatch(fitfunc(formula=formula, data=data,...),
                          error=function(e){NULL})
  original <- list(formula=Reduce(paste, deparse(formula)),
                   criterion=criterion(originalfit),
                   model=originalfit)
 

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
  solutions$criterion <- Cri[[info$solution]]
  solutions$swaps <- info$iteration
  if (isTRUE(return.models)) {
    solutions <- as.tibble(solutions)
    solutions$swapped.to.model <- list(NA)
    solutions$checked.model <- list(NA)
    for (k in 1:numrs) {
      solutions$swapped.to.model[[k]] <- MDL[unique(info$history[[k]])]
      solutions$checked.model[[k]] <- MDL[unique(info$steps[[k]])]
    }
  } else {
    solutions <- data.frame(solutions, stringsAsFactors = FALSE)
    rownames(solutions) <- NULL
  }

  
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
  table <- as_tibble(table)
  if (isTRUE(return.models)) {
    table$model <- MDL[sln.keys]
  } else {
    table$model <- lapply(table$formula, FUN = function(form) {
      fitfunc(form, data=data, ...)
    })
  }


  res <- list(solutions=solutions, table=table,
              nfits=Cri$size(), nchecks=sum(info$check),
              original=original)
  class(res) <- "FSA"
  return(res)
}


#' @export
#' @describeIn FSA alias for \code{FSA(fitfunc=lm,...)}
lmFSA <- function( ...) {FSA(fitfunc = lm, ...)}

#' @export
#' @describeIn FSA alias for \code{FSA(fitfunc=glm,...)}
glmFSA <- function( ...) {FSA(fitfunc = glm, ...)}
