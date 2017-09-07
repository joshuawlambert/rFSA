#' Diagnostic Plots for FSA solutions
#'
#' @param x FSA object to see diagnostic plots on.
#' @param ask logical; if TRUE, the user is asked before each plot. See help(plot.lm). 
#' @param easy logical; should diagnostic plots be presented in easy to read format?
#' @param ... arguments to be passed to other functions. 
#' @return diagnostic plots to plot window. 
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
#' plot(x=fit)
plot.FSA <- function(x,ask = F,easy = T,...) {
  stopifnot(inherits(x, "FSA"))
  fm <- fitmodels(x)
  if (length(fm) < 4) {
    dm <- length(fm)
  } else
    dm <- 4
  if (easy == F) {
    par(mfrow = c(1,4))
  }
  else{
    par(mfrow = c(dm,4))
  }
  for (i in 1:length(fm)) {
    one<-capture.output(pfit<-print(x))
    plot(fm[[i]],ask = ask,main = rownames(pfit)[i])
  }
}






