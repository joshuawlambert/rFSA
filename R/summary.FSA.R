#' Summary function for FSA solutions
#'
#' @param object FSA object to see summaries on.
#' @param ... arguments to be passed to other functions. 
#' @return list of summarized lm or glm output. 
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
#'            quad=FALSE,m=2,numrs=10,save_solutions=FALSE,cores=1)
#' summary(fit)
summary.FSA <- function(object,...) {
  fm <- fitmodels(object)
  sumresls <- list()
  for (i in 1:length(fm)) {
    sumresls[[i]] <- summary(fm[[i]])
  }
  return(sumresls)
}
