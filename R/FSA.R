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

    ##Generate random starting positions
    starts <- replicate(n=numrs, expr=pos2key(sort(sample(xpos, m, replace = F))))

    if(usehist)
    {
        hist.sln <- hash()
        hist.pos <- as.list(starts)
        names(hist.pos) <- starts
    }
    cur.pos <- starts
    sln <- c()
    while(length(cur.pos)>0)
    {
        ##Find stepping positions
        steps<-unique.array(matrix(unlist(lapply(1:length(cur.pos),FUN = function(x){swaps(cur = key2pos(cur.pos[x]),n = P)})),ncol = m,byrow = T),MARGIN = 1)

        ##Calculate criterion for each position
        form <- function(val)
        {
            as.formula(paste0(yname, "~", paste(allname[val], collapse="+")))
        }
        Cri <- hash()
        mclapply(
            X=1:nrow(steps), mc.cores=cores,
            FUN = function(x)
            {
                key <- pos2key(steps[x,])
                if(!has.key(key, Cri))
                    Cri[[key]] <- criterion(method(formula=form(steps[x,]), data = data,...))
            })
        
        ##Find the best next position for each current position
        tmp <- mclapply(
            X=1:length(cur.pos), mc.cores=cores,
            FUN = function(x)
            {
                step <- swaps(key2pos(cur.pos[x]), P)
                criterions <- Cri[
                    sapply(X=1:ncol(step),
                           FUN = function(x){pos2key(step[,x])})]
                keys(criterions)[which.min(values(criterions))]
            })
        next.pos <- unlist(tmp)
        names(next.pos) <- cur.pos
        if(usehist) stopifnot(length(cur.pos)==length(hist.pos))

        ##Check if any solutions are found
        mask <- cur.pos == next.pos
        cur.sln <- rep(NA, length(cur.pos))
        if(any(mask)) cur.sln[mask] <- cur.pos[mask]
        if(usehist)
        {
            mask2 <- has.key(next.pos, hist.sln)
            if(any(mask2))
                cur.sln[mask2] <- values(hist.sln, next.pos[mask2])
            mask <- mask | mask2
        }
        names(cur.sln) <- cur.pos
        sln <- c(sln, cur.sln[mask])


        if(usehist)
        {
            ##Update solution history
            stopifnot(is.vector(cur.sln))
            for(key in cur.pos[mask])
                hist.sln[hist.pos[[key]]] <- cur.sln[key]
        
            ##Update position history
            for(key in cur.pos[!mask])
                hist.pos[[key]] <- c(hist.pos[[key]], next.pos[key])
        }

        ##Update settings and iterate
        if(usehist) hist.pos <- hist.pos[!mask]
        cur.pos <- next.pos[!mask]
        if(usehist) names(hist.pos) <- cur.pos
    }
    return(sln)
}
