#' Prediction function for FSA solutions
#'
#' @param object FSA object to conduct predictions on.
#' @param ... other parameters passed to fitmodels or predict functions. See help(fitmodels) or help(predict) for assistance.
#' @return list of predicted values from each FSA model. 
#' @importFrom graphics par plot
#' @importFrom stats AIC anova as.formula cov fitted formula glm influence lm predict resid var
#' @importFrom utils capture.output tail
#' @export
#'
#' @examples
#' #use mtcars package see help(mtcars)
#' data(mtcars)
#' colnames(mtcars)
#' fit<-lmFSA(formula="mpg~hp*wt",data=mtcars,fixvar="hp",
#'              quad=FALSE,m=2,numrs=10,save_solutions=FALSE,cores=1)
#' predict(fit)
#' predict(fit,newdata=mtcars[1:15,])
predict.FSA <- function(object,...) {
  stopifnot(inherits(object, "FSA"))
  fm <- fitmodels(object,...)
  res <- list()
  for (i in 1:length(fm)) {
    res[[i]] <- predict(fm[[i]],...)
  }
  return(res)
}
