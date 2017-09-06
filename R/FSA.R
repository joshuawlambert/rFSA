#' FSA: A function to find best subsets and interactions in statistical models.
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
#' @importFrom hash hash has.key keys values
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
FSA <- function(formula, data, fitfunc=lm, fixvar = NULL, quad = FALSE,
                m = 2, numrs = 1, cores=1, interactions = T,
                criterion = AIC, minmax="min", checkfeas=NULL, var4int=NULL,
                min.nonmissing=1,...)
{
  formula <- as.formula(formula)
  data <- data.frame(data)

  if(.Platform$OS.type != "unix") cores = 1

  if(!(is.atomic(min.nonmissing) & length(min.nonmissing)==1))
    stop("min.nonmissing should be a scalar.")
  
  yname <- all.vars(formula)[1]
  allname <- colnames(data)
  stopifnot(all(c(yname,fixvar) %in% allname))
  P <- length(allname)-1
  ypos <- which(allname == yname)
  xpos <- setdiff(1:(P+1), ypos)
  criterion <- criterion
  method <- fitfunc
  which.best <- switch(tolower(minmax), min=which.min, max=which.max)
  bad.cri <- switch(tolower(minmax), min=Inf, max=-Inf)

  ##Generate random starting positions
  #if checkfeas != NULL and length(checkfeas)==m then put the the check feas in the last position of starts
  if (is.null(checkfeas))
    starts <- replicate(n=numrs, expr=pos2key(sort(sample(xpos, m, replace = F))))
  else if (length(checkfeas)!=m)
    stop("sorry, the number of variables in checkfeas is not equal to m. Please try again.")
  else
  {
    starts <- replicate(n=numrs, expr=pos2key(sort(sample(xpos, m, replace = F))))
    starts[length(starts)]<-pos2key(which(colnames(data) %in% checkfeas))
  }
  
  ## cur.key stores the keys of current positions
  ## during optimization.
  cur.key <- starts
  sln <- c()
  
  form <- if(interactions==F){function(val){as.formula(paste0(yname, "~", paste(fixvar,collapse = "+"),"+",paste(allname[val], collapse="+")))}
  } else{function(val){as.formula(paste0(yname, "~", paste(fixvar,collapse = "+"),"+", paste(allname[val], collapse="*")))}}

  ## Initilize a hash table to store criterion for the computed
  ## combinations. The keys used to index criterions are
  ## produced by pos2key, and could be decoded by key2pos
  Cri <- hash()
  
  while(length(cur.key)>0)
  {
    ##Find stepping positions
    ## steps<-unique.array(matrix(unlist(lapply(
    ##   1:length(cur.key),
    ##   FUN = function(x){
    ##     if(is.null(var4int)){swaps(cur = key2pos(cur.key[x]), n = P + 1, quad = quad,yindex = ypos)
    ##     } else{swaps(cur = key2pos(cur.key[x]), n = P + 1, quad = quad,yindex = ypos)[,which(apply(swaps(cur = key2pos(cur.key[x]), n = P + 1, quad = quad,yindex = ypos)==which(colnames(data)==var4int),MARGIN = 2,FUN = any))]
    ##     }
    ##   })),ncol = m,byrow = T),MARGIN = 1)
    steps <- lapply(
      cur.key,
      FUN = function(x)
      {
        tmp <- swaps(cur = key2pos(x), n=P+1, quad = quad,yindex = ypos)
        if (!is.null(var4int))
          tmp <- tmp[,which(apply(tmp==which(allname==var4int), MARGIN=2, FUN=any))]
        apply(tmp, MARGIN=2, FUN=pos2key)
      }
    )
    names(steps) <- cur.key

    ## Calculate criterion for each next step position
    ## Basically, we will check the criterion hash table
    ## If a combination's criterion already exists in
    ## the hash table, we simply use it. Otherwise, we
    ## will calculate the criterion and insert it into
    ## the hash table
    steps.noCri <- unique(unlist(steps))
    steps.noCri <- steps.noCri[!has.key(steps.noCri, Cri)];
    if (length(steps.noCri) > 0 )
    {
      new.Cri <- unlist(mclapply(
        steps.noCri, mc.cores=cores,
        FUN = function(key)
        {
          pos <- key2pos(key)
          ## check if there are too many NAs
          if (sum(!apply(is.na(data[,c(pos, ypos)]), 1, any)) < min.nonmissing)
            newCri = bad.cri
          else
            newCri <- tryCatch(criterion(method(formula=form(pos), data = data,...)),
                               error=function(cond) {bad.cri})
          names(newCri) <- key
          return(newCri)
        }))
      Cri[names(new.Cri)] <- new.Cri
    }
    
    ##Find the best next position for each current position
    next.key <- unlist(mclapply(
      cur.key, mc.cores=cores,
      FUN = function(key)
      {
        step <- steps[[key]]
        criterions <- Cri[step]
        return(keys(criterions)[which.best(values(criterions))])
      }
    ))
    names(next.key) <- cur.key

    ##Check if any solutions are found
    mask <- cur.key == next.key
    cur.sln <- rep(NA, length(cur.key))
    cur.sln[mask] <- cur.key[mask]
    names(cur.sln) <- cur.key
    sln <- c(sln, cur.sln[mask])

    ##Update settings and iterate
    cur.key <- next.key[!mask]
  }
  
  #formatting results for output
  for(i in 1:dim(table(sln))){
    prev=paste(lapply(1:dim(table(sln)),FUN=function(x){colnames(data)[eval(parse(text = paste("c(",names(table(sln))[x],")")))]})[[i]],collapse = ", ");
    if(i==1){prev1=NULL;crit=NULL;form.1=NULL}
    prev1=c(prev1,prev)
    form.1<-c(form.1,format(form(eval(parse(text = paste("c(",names(table(sln))[i],")"))))))
    crit<-c(crit,criterion(method(formula=form(eval(parse(text = paste("c(",names(table(sln))[i],")")))), data = data,...)))
  }
  sln_return<-data.frame(starts=starts,solutions=sln)
  rownames(sln_return)<-NULL
  solutions<-data.frame(FS.Num=1:dim(table(sln)),
                        Formula=as.character(form.1),
                        Positions=names(table(sln)),"Variable Names"=prev1,
                        Criterion=crit,
                        "Times"=as.numeric(table(sln))
                        )
  return(list(results=sln_return,solutions=solutions))
}
