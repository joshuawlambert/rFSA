#' Fitted Values for FSA solutions
#'
#' @param object FSA object to get fitted values from.
#' @param ... other parameters passed to fitmodels or fitted function. See help(fitmodels) or help(fitted) for assistance.
#' @return list of fitted values from each FSA model. 
#' @importFrom graphics par plot
#' @importFrom methods show
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
#' fitted(fit)
fitted.FSA <- function(object,...) {
  stopifnot(inherits(object, "FSA"))
  fm <- fitmodels(object,...)
  res <- list()
  for (i in 1:length(fm)) {
    res[[i]] <- fitted(fm[[i]],...)
  }
  return(res)
}