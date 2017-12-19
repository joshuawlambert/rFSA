# Overview
This R package utilizes an exchange algorithm to find best higher order interactions, and best subsets in statistical models. The algorithm searches the data space for models of a specified form that are statistically optimal. Many replications of this algorithm will produce a set of `feasible solutions', which the researcher can investigate. The algorithm can help improve existing models used in bioinformatics, health care, or other fields which have yet to explore quadratic terms, interactions, or a higher order of predictors because of the size of their datasets. The package, rFSA, is flexible to many different statistical methods and criteria functions. Statistical methods in R that have a formula, and data command can be used by rFSA to search for interactions or best subsets. Users can use common criterion functions (R squared, AIC, PRESS, ... ) or write their own [See Guide](http://www.shinyfsa.org)

# The FSA Algorithm
Statisticians are often faced with the problem of identifying a subset of k explanatory variables from p variables Xp, including interactions and quadratic terms. Consider fixing p+ explanatory variables in a preliminary model. Denote these variables Xp+. Let m(Y;Xp+) be an objective function that can be a measure of model quality i.e., R2; AIC; BIC; etc. We
wish to find the k additional variables denoted Xk to add to the model that optimizes the objective function m(Y;Xp+;Xk).

The Feasible Solution Algorithm (FSA) addresses this problem in the following way:
1. Choose Xk randomly and compute the objective function m.
2. Consider exchanging one of the k selected variables from the current model
3. Make the single exchange that improves the objective function m the most.
4. Keep making exchanges until the objective function does not improve. These variables Xp+;Xk are called a feasible solution.
5. Return to (1) to find another feasible solution.

In another instance of the FSA, we include the jth order interaction and lower order terms we are considering in step 1. We then continue on to step 2, only this time when we make an exchange it changes the jth order interaction and the lower interactions and main effects as well. We could then optimize based on a model criterion or on an interaction terms p-value.

A single iteration of FSA yields a feasible solution in the sense that it may globally optimize m(Y;X). Of course, the algorithm may converge somewhere other than the global optimum. Using the algorithm multiple times identifies multiple feasible solutions, the best of which may be the global optimum.

Miller (1984) outlines the FSA described above in the following way. Suppose we have 26 predictor variables labeled A through Z. Imagine you wish to find the best subset with four predictor variables. First start randomly with four predictors. Suppose these are ABCD. Consider changing one of A, B, C, or D with one of the other 22 remaining variables. Make the change that improves the objective function the most. Suppose we swap C for X. Now we have ABXD. Next consider changing one of A, B, X, or D (Considering X here is redundant and not necessary). This process is repeated until no further improvements, to the objective function, can be made.

This method, coupled with repeating it for different random starts, can give different solutions which could be interesting from a clinical or scientific viewpoint. These unique feasible solutions are optimal for the criterion function that was chosen by the user in that no one exchange of any one variable can improve the criterion function.

# The R package rFSA
To use the R package, rFSA, the user must know how they wish to model the data (method) and how the will evaluate the fit of the models that will be checked (criterion function). 
## Example (mtcars)
In R, a commond test dataset to analyze is `mtcars`. For this example we will use multiple linear regression to fit the model and Adjusted R Squared to asses the model fit on the response(Miles Per Gallon).
```R
data(mtcars)
help(mtcars)
```
For this example, let us assume that we have already found that the weight(wt) and the number of cylinders(cyc) in a car are statistically signficant predictors of Miles Per Gallon(mpg) for that car. And, let's say that we wish to explore two-way interactions for inclusion in or model with wt and cyc. 

To do this with rFSA, we could run the following code:
```R
install.packages(rFSA) #or install_github("joshuawlambert/rFSA")
library(rFSA)
data(mtcars)

set.seed(123)

fsaFit<-FSA(
  formula="mpg~wt+cyl", #Model that you wish to compare new models to. The variable to the left of the '~' will be used as the response variable in all model fits
  data=mtcars, #specify dataset 
  fitfunc = lm, #method you wish to use 
  fixvar = c('wt','cyl'), # variables that should be fixed in every model that is considered 
  m = 2, #order of interaction or subset to consider
  numrs = 10, #number of random starts to do 
  interactions = TRUE, #If TRUE, then the m variables under condsideration will be added to the model with a '*' between them, if FALSE then the m variables will be added to the model with a '+' between them. Basically, do you want to look for interactions or best subsets.
  criterion = adj.r.squared, #Criterion function used to asses model fit
  minmax = "max" #Should Criterion function be minimized ('min') or maximized ('max').
  
)

fsaFit #shows results from running FSA
print(fsaFit) #shows results from running FSA
summary(fsaFit) #shows summary from all models found by FSA
plot(fsaFit) #plots diagnostic plots for all models found by FSA
fitted(fsaFit) #fitted values from all models found by FSA
predict(fsaFit) #predicted values from all models found by FSA, can also add newdata command.
```
As we can see from `print(fsaFit)`, from the 10 random starts there were 2 feasible solutions (FS). The two feasible included an interaction between hp*wt and drat*carb.  Each of these FS happened 9 and 1 respectively. After looking at `summary(fsaFit)` we can see that hp*wt is statistically significant (p-value<0.01) and drat*carb is marginally significant (p-value ~= 0.06). If we wished to find interactions that were significant, we could change `criterion=int.p.val` and `minmax = "min"`. Do so, will yeild one FS: hp*wt.

Following up these results with sufficient checks into model fit and diagnositic plots is . 


Other options for the FSA function in rFSA include: 
```R
cores = 1  #*FOR LINUX USERS ONLY* uses parallel package to use multiple cores if available. 
quad = FALSE #should quadratic terms be considered by rFSA (example: wt^2)
checkfeas = NULL #Choose a starting place for FSA. If used, should be a vector same length as m from above. Example: c('wt','cyc')
var4int = NULL #Variable to fix in interaction. Useful when considering 3 or more way interactions.
min.nonmissing = 1 #Don't consider models that have less than or equal to this number of observations
return.models = FALSE #should all models that are checked be returned? Useful when you want to ploc criterion history.
```
see `help(rFSA)' for more details.
