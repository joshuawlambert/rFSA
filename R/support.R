#' Variables to include in first steip of an mth order interaction model determined from the Feasible Soution Alorithm.
#'
#' Finds the swaps available given a current position.
#' @param cur A vector of length greater than 2 of what current explantory varialbes are being used in the model.
#' @param n The number of explanatory variables in available to swap.
#' @param quad Whether to include quadratic terms. ie (x1*x1) as potential swaps.
#' @param yindex index of response variable.
#' @return a matrix with the possible forms by column.
#' @importFrom graphics par plot
#' @importFrom stats AIC anova as.formula cov fitted formula glm influence lm predict resid var
#' @importFrom utils capture.output tail
#' @export
swaps<-function(cur,n,quad=FALSE,yindex){
  m <- length(cur)
  if (!quad) {
    l <- (n - m) * m
    possible <- matrix(rep(cur, l), nrow = l, byrow = T)
    s <- 1:n
    lapply(1:m, FUN = function(i) s <<- s[-which(s == cur[i])])
    lapply(1:m, FUN = function(j) possible[((j - 1) * (n - 
                                                         m) + 1):((j - 1) * (n - m) + (n - m)), j] <<- s)
  }
  if (quad) {
    l <- (n - 1) * m
    possible <- matrix(rep(cur, l), nrow = l, byrow = T)
    s <- 1:n
    smat <- matrix(rep(0, (n - 1) * m), ncol = m)
    lapply(1:m, FUN = function(i) smat[, i] <<- s[-which(s == 
                                                           cur[i])])
    lapply(1:m, FUN = function(j) possible[((j - 1) * (n - 
                                                         1) + 1):((j - 1) * (n - 1) + (n - 1)), j] <<- smat[, 
                                                                                                            j])
  }
  vec <- cbind(cur, t(possible))
  vec<-unique(apply(vec, sort, MARGIN = 2), MARGIN = 2)
  if(length(dim(vec))==0){vec<-t(vec)}
  vec<-vec[,-which(apply(apply(X = vec,MARGIN = 1,function(x){x==yindex}),MARGIN = 1,FUN = any))]
  if(length(dim(vec))==0){vec<-t(vec)}
  return(vec)
}
#' Variables to include in the >1st step of an mth order interaction model determined from the Feasible Soution Alorithm.
#'
#' Finds the swaps available given a current position given previous picks.
#' @param curpos A vector of length greater than 2 of what current explantory varialbes are being used in the model.
#' @param n The number of explanatory variables in available to swap.
#' @param quad Whether to include quadratic terms. ie (x1*x1) as potential swaps.
#' @param prevpos A vector of previous best spots
#' @return a matrix with the possible forms by column.
#' @export
nextswap <- function(curpos,n,prevpos,quad) {
  swps <- swaps(curpos,n,quad)
  nextpos <- rep(FALSE,dim(swps)[2])
  for (i in 1:dim(swps)[1]) {
    nextpos <- nextpos + (swps[i,] %in% prevpos)
  }
  
  retSwps <- swps[,nextpos == (length(curpos) - 2)]
  if (is.null(dim(retSwps))) {
    retSwps <- t(t(retSwps))
  }
  return(list(nswaps = retSwps,prevpos = prevpos))
}

pos2key <- function(pos) { paste(pos,collapse = ",")}
key2pos <- function(key) { eval(parse(text=paste0("c(", key, ")")))}
