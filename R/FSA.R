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
#' @param usehist use history to potentially save computational time.
#' @param ... other arguments passed to fitfunc.
#'
#' @importFrom hash hash
#' @importFrom parallel mclapply
#' @return matrix of results
#' @export
#'
#' @examples
#' 
#' N <- 100 #number of obs
#' P <- 100 #number of variables
#' data <- data.frame(matrix(rnorm(N*(P+1)), nrow = N, ncol = P+1))
#' 
#' sln <- FSA(formula = "X101~1", data = data, cores = 1, m = 2,
#' interactions = F, criterion = AIC, minmax = "min",
#' numrs = 10,usehist=F)
#' sln
FSA <- function(formula, data, fitfunc=lm, fixvar = NULL, quad = FALSE, m = 2,
                numrs = 1, cores=1, interactions = T, criterion = AIC,
                minmax="min", usehist = T,...)
{
    formula <- as.formula(formula)
    data <- data.frame(data)

    if(.Platform$OS.type != "unix") cores = 1

    yname <- all.vars(formula)[1]
    allname <- colnames(data)
    stopifnot(all(c(yname,fixvar) %in% colnames(data)))
    P <- length(allname)-1
    ypos <- which(allname == yname)
    xpos <- setdiff(1:(P+1), ypos)
    criterion <- criterion
    method <- fitfunc
    which.best <- switch(tolower(minmax), min=which.min, max=which.max)

    ##Generate random starting positions
    starts <- replicate(n=numrs, expr=pos2key(sort(sample(xpos, m, replace = F))))

    if(usehist)
    {
        hist.sln <- hash()
        hist.key <- as.list(starts)
        names(hist.key) <- starts
    }
    cur.key <- starts
    sln <- c()
    
    form <- if(interactions==F){function(val){as.formula(paste0(yname, "~", paste(fixvar,collapse = "+"),"+",paste(allname[val], collapse="+")))}
      } else{function(val){as.formula(paste0(yname, "~", paste(fixvar,collapse = "+"),"+", paste(allname[val], collapse="*")))}}

    Cri <- hash() # Initilize criterion records
    
    while(length(cur.key)>0)
    {
        ##Find stepping positions
        steps<-unique.array(matrix(unlist(lapply(1:length(cur.key),FUN = function(x){swaps(cur = key2pos(cur.key[x]),n = P+1,yindex=ypos)})),ncol = m,byrow = T),MARGIN = 1)

        ##Calculate criterion for each position

        ## Modified by Liyu Gong 7/6/2017
        ## Initilize of criterion records should be put
        ## out of the while loop
        ##Cri <- hash()

        ##************************************************************
        ## Modified by Liyu Gong @ 7/6/2017
        ## This implementation to fix the bug is good, but we still
        ## need to deal with the situation when tmp is NULL. So I
        ## modify it
        ##------------------------------------------------------------
        ## tmp<-data.frame(matrix(unlist(mclapply(
        ##     X=1:nrow(steps), mc.cores=cores,
        ##     FUN = function(x)
        ##     {
        ##         key <- pos2key(steps[x,])
        ##         if(!has.key(key, Cri))
        ##             c(criterion(method(formula=form(steps[x,]), data = data,...)),key)
        ##     })),ncol=2,byrow=T))
        ## tmp$X2<-as.character(tmp$X2); tmp$X1<-as.numeric(as.character(tmp$X1));
        ## tmp<-lapply(1:dim(tmp)[1],FUN = function(x){Cri[[tmp$X2[x]]]<-tmp$X1[x]})
        tmp <- mclapply(
            X=1:nrow(steps), mc.cores=cores,
            FUN = function(x)
            {
                key <- pos2key(steps[x,])
                if(!has.key(key, Cri))
                {
                    newCri <- criterion(method(formula=form(steps[x,]), data = data))
                    names(newCri) <- key
                    return(newCri)
                }
                else return(NULL)
            })
        tmp <- unlist(tmp)
        if(!is.null(tmp))
            Cri[names(tmp)] <- tmp
        
        ##Find the best next position for each current position
        tmp <- mclapply(
            X=1:length(cur.key), mc.cores=cores,
            FUN = function(x)
            {
                step <- swaps(key2pos(cur.key[x]), P+1,yindex=ypos)
                criterions <- Cri[
                    sapply(X=1:ncol(step),
                           FUN = function(x){pos2key(step[,x])})]
                keys(criterions)[which.best(values(criterions))]
            })
        next.key <- unlist(tmp)
        names(next.key) <- cur.key
        if(usehist) stopifnot(length(cur.key)==length(hist.key))

        ##Check if any solutions are found
        mask <- cur.key == next.key
        cur.sln <- rep(NA, length(cur.key))
        if(any(mask)) cur.sln[mask] <- cur.key[mask]
        if(usehist)
        {
            mask2 <- has.key(next.key, hist.sln)
            if(any(mask2))
                cur.sln[mask2] <- values(hist.sln, next.key[mask2])
            mask <- mask | mask2
        }
        names(cur.sln) <- cur.key
        sln <- c(sln, cur.sln[mask])

        if(usehist)
        {
            ##Update solution history
            stopifnot(is.vector(cur.sln))
            for(key in cur.key[mask])
                hist.sln[hist.key[[key]]] <- cur.sln[key]
        
            ##Update position history
            for(key in cur.key[!mask])
                hist.key[[key]] <- c(hist.key[[key]], next.key[key])
        }

        ##Update settings and iterate
        if(usehist) hist.key <- hist.key[!mask]
        cur.key <- next.key[!mask]
        if(usehist) names(hist.key) <- cur.key
    }
    
    #formatting results for output
    for(i in 1:dim(table(sln))){
      prev=paste(lapply(1:dim(table(sln)),FUN=function(x){colnames(data)[eval(parse(text = paste("c(",names(table(sln))[x],")")))]})[[i]],collapse = ", ");
      if(i==1){prev1=NULL;crit=NULL;form.1=NULL}
      prev1=c(prev1,prev)
      form.1<-c(form.1,format(form(eval(parse(text = paste("c(",names(table(sln))[i],")"))))))
      crit<-c(crit,criterion(method(formula=form(eval(parse(text = paste("c(",names(table(sln))[i],")")))), data = data,...)))
    }
    
    solutions<-data.frame("FS Num"=1:dim(table(sln)),
                          Formula=form.1,
                          Positions=names(table(sln)),"Variable Names"=prev1,
                          Criterion=crit,
                          "Times"=as.numeric(table(sln))
                          )
    return(solutions)
}
