#' Model fitting function for FSA solutions
#'
#' @param object  FSA object to construct models on.
#' @param ... other parameters passed to lm or glm. See help(lm) or help(glm) for other potential arguements.
#' @return list of FSA models that have been fitted.
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
#'      quad=FALSE,m=2,numrs=10,save_solutions=FALSE,cores=1)
#' fitmodels(fit)
fitmodels <- function(object,...) {
  stopifnot(inherits(object, "FSA"))
  resls <- c(list(object$original$model), object$table$model)
  ## one<-capture.output(forms <- print(object))
  ## if (is.null(object$call$fam)) {
  ##   for (i in 1:(dim(object$table)[1] + 1)) {
  ##     resls[[i]] <- lm(forms$Formula[[i]],data = object$call$data,...)
  ##   }
    
  ## } else{
  ##   for (i in 1:(dim(object$table)[1] + 1)) {
  ##     resls[[i]] <-
  ##       glm(forms$Formula[[i]],data = object$call$data,family = object$call$fam,...)
  ##   }
  ## }
  return(resls)
}
