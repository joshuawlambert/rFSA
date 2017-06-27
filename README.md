# Overview
This R package applies a new algorithm intended to improve statistical models for prediction. The algorithm searches the data space for models of a specified form that are statistically optimal. Many replications of this algorithm will produce a set of `feasible solutions', which the researcher can investigate. The algorithm can help improve existing models used in bioinformatics, health care, or other fields which have yet to explore quadratic terms, interactions, or a higher order of predictors because of the size of their datasets. The package, rFSA, is exible to many different model forms and criteria functions. Currently, linear models and generalized linear models are supported.

#FSA Algorithm Details
Data analysts are often faced with the problem of identifying a subset of k explanatory variables from p variables Xp, including interactions and quadratic terms. Consider fixing p+ explanatory variables in a preliminary model. Denote these variables Xp+. Let m(Y;Xp+) be an objective function that can be a measure of model quality i.e., R2; AIC; BIC; etc. We
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