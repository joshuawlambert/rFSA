#' twFSA
#'
#' @export
#' 
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
