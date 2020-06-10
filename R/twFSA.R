#' twFSA
#'
#' @description A function for termwise feasiblity
#' 
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. 
#' @param data a data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
#' @param fitfunc the method that should be used to fit the model. For Example: lm, glm, or other methods that rely on formula, data, and other inputs.
#' @param fixvar variable(s) to fix in the model. Usually a covariate that should always be included (Example: Age, Sex). Will still consider it with interactions. Default is NULL.
#' @param quad Include quadratic terms or not. Logical.
#' @param cores number of cores to use while running. Note: Windows can only use 1 core. See mclapply for details. If function detects a Windows user it will automatically set cores=1.
#' @param criterion which criterion function to either maximize or minimize. For linear models one can use: r.squared, adj.r.squared, cv5.lmFSA (5 Fold Cross Validation error), cv10.lmFSA (10 Fold Cross Validation error), apress (Allen's Press Statistic), int.p.val (Interaction P-value), AIC, BIC.
#' @param minmax whether to minimize or maximize the criterion function
#' @param checkfeas vector of variables that could be a feasible solution. These variables will be used as the last random start.
#' @param var4int specification of which variables to check for marginal feasiblilty. Default is NULL
#' @param min.nonmissing the combination of predictors will be ignored unless this many of observations are not missing
#' @param ... other arguments passed to fitfunc.
#'
#' @import hash
#' @importFrom parallel mclapply
#' @import tibble
#' @return matrix of results
#' @export
twFSA<-function(formula, data, fitfunc=lm, fixvar = NULL, quad = FALSE,    
                cores=1,
                criterion = AIC, minmax="min", checkfeas=NULL, var4int=NULL,
                min.nonmissing=1, ...){
  
  tmp <- as.character(as.formula(formula))
  yname <- tmp[2]
  terms <- strsplit(tmp[3], " [+] ")[[1]]
  
  for (i in 1:length(terms)) {
    formula.1 <- paste0(yname, "~.")
    
    ## construct fixed part
    fix.formula.1 <- paste0(terms[-i], collapse = " + ")
    
    ## deduce m from unfixed term
    m <- length(all.vars(paste0("~", terms[i])))
    
    ## deduce interaction from unfixed term
    interaction <- grepl("*", terms[i])
    
    ## call FSA
    res <- FSA(formula=formula.1, data=data, fitfunc=fitfunc, fixvar=fixvar,
               m = 2, interactions=interaction, fix.formula=fix.formula.1,
               quad=quad, cores=cores, criterion=criterion, minmax=minmax,
               checkfeas=checkfeas, var4int=var4int,
               min.nonmissing=min.nonmissing, return.models = FALSE)
    
    ## extract best term found by FSA and replace the corresponding term
    tmp <- as.character(as.formula(res$table$formula))[[3]]
    terms[i] <- tail(strsplit(tmp, split=' [+] ')[[1]], n=1)
  }
  
  return(res)
}