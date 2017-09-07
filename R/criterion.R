#' @export
globalVariables(c("newdata","moves","solutions","myWarnings"))

#' An rFSA Criterion Function.
#' @description  rFSA Criterion Function to compute R squared.
#' @param model lm or glm fit to be passed. 
#' @param name passed to print.FSA
#' @importFrom graphics par plot
#' @importFrom methods show
#' @importFrom stats AIC anova as.formula cov fitted formula glm influence lm predict resid var
#' @importFrom utils capture.output tail
#' @export
r.squared <- function(model,name = "R Squared") {
  #R squared
  if(is.null(model$family[[1]])){return(summary(model)$r.squared)
    } else return(1.1)
}

#' An rFSA Criterion Function.
#' @description  rFSA Criterion Function to compute Root Mean Squared Error.
#' @param model lm or glm fit to be passed. 
#' @param name passed to print.FSA
#' @importFrom graphics par plot
#' @importFrom methods show
#' @importFrom stats AIC anova as.formula cov fitted formula glm influence lm predict resid var
#' @importFrom utils capture.output tail
#' @export
rmse <- function(model,name = "RMSE") {
  #Root Mean Squared Error
  sqrt(mean(model$residuals ^ 2))
}

#' An rFSA Criterion Function.
#' @description  rFSA Criterion Function to compute Adjusted R-Squared.
#' @param model lm or glm fit to be passed. 
#' @param name passed to print.FSA
#' @importFrom graphics par plot
#' @importFrom methods show
#' @importFrom stats AIC anova as.formula cov fitted formula glm influence lm predict resid var
#' @importFrom utils capture.output tail
#' @export
adj.r.squared <- function(model,name = "Adj R Squared") {
  #R squared
  if(is.null(model$family[[1]])){return(summary(model)$adj.r.squared)
  } else return(1.1)
}

#' An rFSA Criterion Function.
#' @description  rFSA Criterion Function to Allen's Press Statistic.
#' @param model lm or glm fit to be passed. 
#' @param name passed to print.FSA
#' @importFrom graphics par plot
#' @importFrom methods show
#' @importFrom stats AIC anova as.formula cov fitted formula glm influence lm predict resid var
#' @importFrom utils capture.output tail
#' @export
apress <- function(model, name = "PRESS") {
  #Allen's PRESS statistic
  presid <- resid(model)/(1 - influence(model)$hat)
  sum(presid^2)
}

#' An rFSA Criterion Function.
#' @description  rFSA Criterion Function to compute Liklihood Ratio Test Statistics p-value for the largest order interation term.
#' @param model lm or glm fit to be passed. 
#' @param name passed to print.FSA
#' @importFrom graphics par plot
#' @importFrom methods show
#' @importFrom stats AIC anova as.formula cov fitted formula glm influence lm predict resid var
#' @importFrom utils capture.output tail
#' @export
int.p.val <- function(model,name = "Interaction P-Value") {
  if((is.null(model$call$family)|is.null(model$family[[1]])) & !(length(grep(pattern = ":",names(model$coefficients)))==0)){return(tail(anova(model,test="LRT")[,5],2)[1])
  }  else if((is.null(model$call$family)|is.null(model$family[[1]])) & (length(grep(pattern = ":",names(model$coefficients)))==0)){
      return(0)
  } else if((model$call$family=="binomial"|model$family[[1]]=="binomial") & !(length(grep(pattern = ":",names(model$coefficients)))==0)){
      return(tail(anova(model,test = "LRT")[,5],1))
  } else if(!(length(grep(pattern = ":",names(model$coefficients)))==0)) {
      return(max(summary(model)$coef[grep(pattern = ":",names(model$coefficients)),4]))
    } else return(0)
}

#' An rFSA Criterion Function.
#' @description  rFSA Criterion Function to compute the Bhattacharyya distance.
#' @param model lm or glm fit to be passed. 
#' @param name passed to print.FSA
#' @importFrom graphics par plot
#' @importFrom methods show
#' @importFrom stats AIC anova as.formula cov fitted formula glm influence lm predict resid var
#' @importFrom utils capture.output tail
#' @export
bdist <- function(model,name = "B Distance") {
  if(length(grep(":",names(model$coefficients)))==0){return(0)} #no interaction in model to check
  if(is.null(model$family)){return(0)} #not logistic regression
  
  nam<-c(all.vars(model$formula)[1],strsplit(names(model$coefficients)[grep(":",names(model$coefficients))],"[:]")[[1]])
  tmp_dat <- eval(model$model)
  y <- tmp_dat[,nam[1]]
  x1 <- tmp_dat[,nam[2]]
  x2 <- tmp_dat[,nam[3]]
  
  X <- cbind(x1,x2)
  a <- X[which(y == 0),]
  b <- X[which(y == 1),]
  
  mu11 <- mean(a[,1])
  mu12 <- mean(a[,2])
  mu21 <- mean(b[,1])
  mu22 <- mean(b[,2])
  
  var1_11 <- var(a[,1])
  var1_22 <- var(a[,2])
  var1_12 <- cov(a[,1],a[,2])
  
  var2_11 <- var(b[,1])
  var2_22 <- var(b[,2])
  var2_12 <- cov(b[,1],b[,2])
  
  mu1 <- c(mu11,mu12)
  mu2 <- c(mu21,mu22)
  sig1 <- matrix(c(var1_11,var1_12,var1_12,var1_22),ncol = 2)
  sig2 <- matrix(c(var2_11,var2_12,var2_12,var2_22),ncol = 2)
  
  sig <- (sig1 + sig2) / 2
  a <- sig[1,1]
  b <- sig[1,2]
  c <- sig[2,1]
  d <- sig[2,2]
  temp <- matrix(c(d,-c,-b,a),ncol = 2)
  sig.inv <- 1 / (a * d - b * c) * temp
  distance <-
    1 / 8 * (t(mu1 - mu2) %*% sig.inv %*% (mu1 - mu2)) + 1 / 2 * log(det(sig) /
                                                                       sqrt(det(sig1) * det(sig2)))
  return(c(distance))
}

#' An rFSA Internal Function.
#' @description  rFSA function to compute the maximum value from a vector with NA's. 
#' @param vec Vector to be passed.
#' @importFrom graphics par plot
#' @importFrom methods show
#' @importFrom stats AIC anova as.formula cov fitted formula glm influence lm predict resid var
#' @importFrom utils capture.output tail 
#' @export
which.max.na <- function(vec) {
  maxval <- max(vec,na.rm = T)
  which(vec == maxval)
}

#' An rFSA Internal Function.
#' @description  rFSA function to compute the minimum value from a vector with NA's. 
#' @param vec Vector to be passed. 
#' @importFrom graphics par plot
#' @importFrom methods show
#' @importFrom stats AIC anova as.formula cov fitted formula glm influence lm predict resid var
#' @importFrom utils capture.output tail
#' @export
which.min.na <- function(vec) {
  minval <- min(vec,na.rm = T)
  which(vec == minval)
}

#' List all included Criteria function for lmFSA and glmFSA.
#'
#' @return list of functions and whether lmFSA or glmFSA work with those functions.
#' @importFrom graphics par plot
#' @importFrom methods show
#' @importFrom stats AIC anova as.formula cov fitted formula glm influence lm predict resid var
#' @importFrom utils capture.output tail 
#' @export
#'
#' @examples
#' list.criterion()
list.criterion<-function(){
  show("Accepted Criteria Functions for lmFSA and glmFSA")
  show("")
  tab<-rbind(c("r.squared","lmFSA"),c("adj.r.squared","lmFSA"),c("AIC","lmFSA; glmFSA"),c("BIC","lmFSA; glmFSA"),
        c("rmse","lmFSA; glmFSA"),c("apress","lmFSA; glmFSA"),c("int.p.val","lmFSA; glmFSA"),
        c("bdist","glmFSA (only 2 way interactions)"))
  colnames(tab)<-c('Criterion',"Accepted in")
  tab<-data.frame(tab)
  show(tab)
  show("")
  show("You can write your own criterion function too! Or use other criterion functions from other packages. Just follow the standard format used in int.p.val or apress as an example.")
}


