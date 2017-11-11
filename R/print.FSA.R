#' Printing function for FSA solutions
#'
#' @param x FSA object to print details about. 
#' @param ... arguments to be passed to other functions. 
#' @return list of Feasible Solution Formula's, Original Fitted model formula and criterion function and times converged to details. 
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
#'                    quad=FALSE,m=2,numrs=10,save_solutions=FALSE,cores=1)
#' print(fit)
print.FSA <- function(x,...)
{
  stopifnot(inherits(x,"FSA"))
  tab <- x$table
  tab$model <- NULL
  tab$opt.criterion <- lapply(tab$opt.criterion, FUN = function(x){paste0(x,collapse=",")})
  tab <- as.data.frame(tab)

  ## drop columns named with Var*
  tab <- tab[,-grep("^Var",names(tab))]

  ## add a row for original fit
  original <- x$original
  original$times <- NA
  original$model <- NULL
  original$criterion <- NA
  original$opt.criterion <- NA
  
  original <- as.data.frame(original, stringsAsFactors = FALSE)
  tab <- rbind(original, tab)
  if (length(grep("^criterion.",names(tab))) == 1) {
##  if (sum(startsWith(names(tab),"criterion."))==1) {
    tab$criterion <- NULL
    tab$opt.criterion <- NULL
    names(tab)[names(tab)=="criterion.1"] <- "criterion"
  }
  rownames(tab) <- c("Original Fit", paste0("FS", 1:(nrow(tab)-1)))
  print(tab)
}

## print3.FSA <- function(x,...) {
##   stopifnot(inherits(x, "FSA"))
##   ## vars <-
##   ##   x$table[,-which(colnames(x$table) %in% c("criterion","times","fixvar"))]
##   vars <- x$table[,-which(colnames(x$table) %in% c("criterion","times","fixvar","formula"))]
##   orgvars <- all.vars(x$call$formula)
##   if (x$call$interactions == T) {
##     form <-
##       function(j)
##         paste0(orgvars[1]," ~ ",if (!is.null(x$call$fixvar)) {
##           paste0(x$call$fixvar,collapse = " + ")
##         },if (!is.null(x$call$fixvar)) {
##           " + "
##         },paste(vars[j,],collapse = " * "),sep = "")
##   }
##   if (x$call$interactions == F) {
##     form <-
##       function(j)
##         paste0(orgvars[1]," ~ ",if (!is.null(x$call$fixvar)) {
##           paste0(x$call$fixvar,collapse = " + ")
##         },if (!is.null(x$call$fixvar)) {
##           " + "
##         },paste(vars[j,],collapse = " + "),sep = "")
##   }
##   tab <- cbind(matrix(lapply(
##     X = 1:dim(x$table)[1],FUN = form
##   )),x$table$criterion,x$table$times)
##   tab <- rbind(c(
##     Reduce(paste0,deparse(formula(x$originalfit))),x$call$criterion(x$originalfit),NA
##   ),tab)
##   tab <- data.frame(tab)
##   cname <- formals(x$call$criterion)$name
##   if (is.null(cname)) {
##     cname = "criterion"
##   }
##   colnames(tab) <- c("Formula", cname, "Times FS")
##   tab[,2] <- as.numeric(unlist(tab[,2]))
##   rownames(tab) <-
##     c("Original Fit   ",paste("FS",1:dim(x$table)[1],sep = ""))
##   tab <- data.frame(tab)
##   print(tab)
## }
