# SurrogateMinimalDepth
In this R package functions are provided to select important variables with surrogate minimal depth (SMD) (main function: var.select.smd) and mimimal depth (MD) (main function: var.select.md) and to investigate variable relations  with the mean adjusted agreement of surrogate variables (main function: var.relations). 

# Install
library(devtools)
install_github("StephanSeifert/SurrogateMinimalDepth")

# Example data
Example data to apply the functions in this package is provided: 

The dataset contains 100 individuals (rows), 1 dependent variable (first column, y), and 200 additional variables (columns 2 to 201): For the simulation of the 200 additional variables nine variables X1,…,X9 called basic variables were simulated. X1 to X6 are causal with constant effect sizes of 1 and X7 to X9 are noncausal with effect size of 0. The outcome y is a linear combination of the causal predictor variables and a normally distributed error term. All basic variables were sampled from a normal distribution (N(0,1)) just like the noise (N(0,0.2)). Furthermore, for three causal and three non-causal variables correlated variables were simulated. Hence, for each of the six basic variables X1, X2, X3, X7, X8, and X9 ten variables with predefined correlations of 0.9 for X1 and X7, 0.6 for X2 and X8, and 0.3 for X3 and X9 were obtained by simulateModule function of the R package WGCNA. The ten variables of each basis variable are labeled: Cp_basicvariable_number. Additional non-correlated and independent predictor variables (cgn) were simulated using the standard normal distribution to reach a total number of 200 variables.

# Usage
## Minimal depth
We want to analyze the example data with minimal depth. 
First the package and the data are loaded:
```
library(SurrogateMinimalDepth)
data("SMD_example_data")
```
Variable selection with var.select.md is conducted:
```
res = var.select.md(x=SMD_example_data[,2:ncol(SMD_example_data)],y=SMD_example_data[,1],ntree=1000)
res$var
 [1] "X2"     "X3"     "X4"     "X5"     "X6"     "cp2_7"  "cp3_4"  "cp3_6"  "cgn_68" "cgn_72"
```

x: Matrix or data.frame of predictor variables with variables in columns and samples in rows. 

y: Vector with values of phenotype variable.

ntree: Number of treesfor random forest.

The selected variables are stored in res$var. In this analysis the causal basic variables X2 to X6, as well as the variables cp2_7, cp3_4, and cp3_6 that are correlated to causal variables, and the non-causal variables cgn_68 and cgn_72 are selected. This means the minimal depth values of these variables are below the threshold.  

## Surrogate Minimal depth
Now we want to analyze the example data with surrogate minimal depth. 
First the package and the data are loaded:
```
library(SurrogateMinimalDepth)
data("SMD_example_data")
```
Variable selection with var.select.smd is conducted:
```
res = var.select.smd(x=SMD_example_data[,2:ncol(SMD_example_data)],y=SMD_example_data[,1],s=10,ntree=1000)
res$var
 [1] "X1"     "X2"     "X3"     "X4"     "X5"     "X6"     "X8"     "cp1_1"  "cp1_2"  "cp1_3"  "cp1_4"  "cp1_5"  "cp1_6"  "cp1_7"  "cp1_8" 
[16] "cp1_9"  "cp1_10" "cp2_1"  "cp2_3"  "cp2_4"  "cp2_6"  "cp2_7"  "cp2_9"  "cp2_10" "cp3_4" 
```
x: Matrix or data.frame of predictor variables with variables in columns and samples in rows. 

y: Vector with values of phenotype variable.

ntree: Number of treesfor random forest.

s: Predefined number of surrogate splits 

The selected variables are stored in res$var. In this analysis the causal basic variables X1 to X6, as well as the variables cp1_1 to cp1_10, cp2_1, cp2_3, cp2_4, cp2_6, cp2_7, cp2_9, cp2_10, and cp3_4" that are correlated to causal variables are selected. This means the minimal depth values of these variables are below the threshold. 

## Variable relations (based on the mean adjusted agreement of surrogate variables)
Now we want to analyze the relations between the variables in the example data. We want to know which variables of the first 100 variables are related to X1 and X7. (Quick reminder: The data set contains 10 correlated variables for each of the two basic variables.) . 
One possibility to investigate variable relations is to used the results from var.select.smd. Hence, first SMD is conducted like in the previous section and the names of the variables are extracted from the result:

```
library(SurrogateMinimalDepth)
data("SMD_example_data")
res = var.select.smd(x=SMD_example_data[,2:ncol(SMD_example_data)],y=SMD_example_data[,1],s=10,ntree=1000)
allvariables=names(res$info$depth)
```
Subsequently, variable relations are analyzed with var.relations:

```
rel=var.relations(trees=res$trees,variables=c("X1","X7"),allvariables=allvariables,candidates=allvariables[1:100],t=5)
rel$var
$X1
 [1] "cp1_1"  "cp1_2"  "cp1_3"  "cp1_4"  "cp1_5"  "cp1_6"  "cp1_7"  "cp1_8"  "cp1_9"  "cp1_10"

$X7
 [1] "cp7_1"  "cp7_2"  "cp7_3"  "cp7_4"  "cp7_5"  "cp7_6"  "cp7_7"  "cp7_8"  "cp7_9"  "cp7_10"
```
trees:	List of trees that was generated by var.select.smd function. 

variables:	Variable names for which related variables should be searched for (has to be contained in allvariables).

allvariables: Vector of all variable names in the original data set.

candidates: Vector of variable names that are candidates to be related to the variables (has to be contained in allvariables).

t: Variable to calculate threshold.

All of the variables that are correlated to X1 are correcly identified as related to X1 and all of the variables that are correlated to X7 are correcly identified as related to X7. 
