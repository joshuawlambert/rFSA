[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/rFSA)](https://cran.r-project.org/package=rFSA)
[![Downloads](http://cranlogs.r-pkg.org/badges/rFSA)](http://cranlogs.r-pkg.org/badges/rFSA)
[![Total Downloads](http://cranlogs.r-pkg.org/badges/grand-total/rFSA)](http://cranlogs.r-pkg.org/badges/grand-total/rFSA)

# Overview
Our article ["rFSA: An R Package for Finding Best Subsets and Interactions"](https://journal.r-project.org/archive/2018/RJ-2018-059/index.html) was published in the R Journal in December 2018. The article outlines The Feasible Solution Algorithm (FSA) and the R package rFSA.

# rFSA
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
install.packages(rFSA) #or devtools::install_github("joshuawlambert/rFSA")
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
As we can see from `print(fsaFit)`, from the 10 random starts there were 2 feasible solutions (FS). The two feasible included an interaction between hp*wt and drat*carb.  Each of these FS happened 9 and 1 respectively. After looking at `summary(fsaFit)` we can see that hp and wt interaction is statistically significant (p-value<0.01) and drat and carb interaction is marginally significant (p-value ~= 0.06). If we wished to find interactions that were significant, we could change `criterion=int.p.val` and `minmax = "min"`. Doing so, will yeild one FS: hp*wt.

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
see `help(FSA)' for more details.

## Visualizing Interactions
Visualizing interactions can be quite difficult depending on the types of variables that are involved in the relationship. The goal of *rFSA* is not to assist the user in visualizing the interaction, but the authors recognize that visual tools are often quite useful conveying the results from a statistical model. 

We have found the *sjPlot* package to be very useful for plotting 2-way interactions. More information about the *sjPlot* package can be found here: http://www.strengejacke.de/sjPlot/. 

For the mtcars example above, the two interactions that were found can be plotted very easy using the *sjPlot* function *plot_model* with the *type="int"* option. Below is some example code:
```R
library(sjPlot)
library(rFSA)
fit<-rFSA::fitmodels(fsaFit)
sjPlot::plot_model(fit[[2]],type = "int")
sjPlot::plot_model(fit[[3]],type = "int")
```
![Feasible solution 1](https://github.com/joshuawlambert/Presentations/blob/master/FS1.png)
![Feasible solution 2](https://github.com/joshuawlambert/Presentations/blob/master/FS2.png)
