#' pFSA: Pareto Feasible Solution Algorithm
#'
#' @description A function using a Feasible Solution Algorithm to estimate a set of models which are on the Pareto frontiers for chosen criteria
#' 
#' 
#' @param numFronts integer number of estimated frontiers to return
#' @param pselExpr expression used by function psel to estimate pareto frontiers. help(psel).
#' @param plot.it TRUE/FALSE for whether to plot the pareto frontiers
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
#' @param fix.formula ...
#' @param ... see arguments taken by function FSA or other functions. help(FSA).
#'
#' @import hash
#' @importFrom parallel mclapply
#' @importFrom graphics legend
#' @import tibble
#' @import rPref
#' @import tidyr
#' @return list of a matrix of all models obtained from FSA (fits) and their criteria. Also a matrix of the estimated frontiers that were requested. The Key column in fits, and pbound refers to the column number of the variables contined in the model fit. For instance, Key="42,96" would refer to the model which contains the variable in the 42nd column and 96th column of the designated dataset.
#' 
#' @export
#'
#' @examples
#'\donttest{
#'N <- 1000 #number of obs
#'P <- 100 #number of variables
#'data <- data.frame(matrix(rnorm(N*(P+1)), nrow = N, ncol = P+1))
#'sln <- pFSA(formula = "X101~1", data = data, m = 2,  criterion = c(max_abs_resid,r.squared),
#'    minmax = c("min","max"),numrs = 10,numFronts = 2,
#'    pselExpr =rPref::low(max_abs_resid)*rPref::high(r.squared),plot.it = TRUE)
#'    }
pFSA <- function(numFronts=2,pselExpr=NULL,plot.it=TRUE,formula, data, fitfunc=lm, fixvar = NULL, quad = FALSE,
                 m = 2, numrs = 1, cores=1, interactions = T,
                 criterion = AIC, minmax="min", checkfeas=NULL, var4int=NULL,
                 min.nonmissing=1, return.models=FALSE, fix.formula=NULL,...)
{ 
  if (length(criterion)<2) {
    stop("for Pareto Optimality you need atleast two criteria functions")
  }
  k<- NULL 
  fsaFit<-FSA(formula, data, fitfunc=fitfunc, fixvar=fixvar, quad=quad,
              m=m, numrs=numrs, cores=cores, interactions=interactions,
              criterion=criterion, minmax=minmax, checkfeas=checkfeas,
              var4int=var4int, min.nonmissing=min.nonmissing,
              return.models=return.models, fix.formula=fix.formula,...)
  
  fits<-spread(fsaFit$criData,key = "k",value = "Values")
  fits2<-fits

  ans<-matrix(data = unlist(mclapply(X = 1:dim(fits2)[1],mc.cores = cores,FUN = function(x){
    if(interactions==TRUE){int="*"} else {int="+"}
    form<-as.formula(paste(all.vars(as.formula(formula))[1],"~",paste(colnames(data)[eval(parse(text=paste0("c(", fits[x,1], ")")))],collapse= int)))
    fit_tmp<-fitfunc(formula=form, data=data,...)
    tmp<-NULL
    for(i in 1:(length(criterion))){
      tmp<-c(tmp,criterion[[i]](fit_tmp))
    }
    tmp
  })),byrow = TRUE,ncol = 2)
  fits2[,-1]<-ans
  
  fsaFit<-FSA(formula, data, fitfunc=lm, fixvar=NULL, quad=FALSE,
              m=m, numrs=numrs, cores=1, interactions=FALSE,
              criterion=criterion, minmax=minmax, checkfeas=NULL,
              var4int=NULL,
              return.models=FALSE, fix.formula=NULL)
  
  fits<-spread(fsaFit$criData,key ="k",value = "Values")
  fits2<<-fits
  
  l<-mclapply(X = which(apply(X = fits2, MARGIN = 1, function(x){any(is.na(x))})),
              mc.cores = cores,
              FUN = function(x,...){
                if(interactions==TRUE){int="*"} else {int="+"}
                form<-as.formula(paste(all.vars(as.formula(formula))[1],"~",
                                       paste(colnames(data)[eval(parse(text=paste0("c(", fits2[x,1], ")")))],
                                             collapse= int))
                )
                fit_tmp<-fitfunc(formula=form, data=data)
                for(i in 1:(length(criterion))){
                  fits2[x,-1][,i]<<-criterion[[i]](fit_tmp)
                }
              }
  )
  
  fits3<-fits2
  
  cname<-gsub(pattern = "high|low|[(]|[])]| ",replacement = "",x = as.character(pselExpr))
  cname<-unlist(strsplit(cname,split = "[*]"))
  
  colnames(fits3)<-c("Key",cname)
  pbound<-psel(df = fits3,pselExpr,top_level=numFronts)
  
  if(length(criterion)>2 & plot.it==TRUE){
    "Sorry, plots cannot be made for more than 2 criteria."
  } else{if(plot.it==TRUE){
    par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
    plot(x = pbound[,2],y = pbound[,3],col=pbound$.level,xlab=cname[1],ylab=cname[2],pch=20,
         main = "Estimated Pareto Frontier Graph for \nChosen Criteria and Number of Fronts")
    legend("bottomright", legend=c("Front 1","Front 2","Not Pareto"), col=c(1,2,3),pch=20, title="Est Pareto Front")
  }
  }
  
  return(list(fits=fits2,pbound=pbound))
}
