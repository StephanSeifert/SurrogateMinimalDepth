# SurrogateMinimalDepth
In this R package functions are provided to select important variables with surrogate minimal depth (SMD) and minimal depth (MD) and to investigate variable relations with the mean adjusted agreement of surrogate variables. 

Please cite the following manuscript if you use the package:
S. Seifert, S. Gundlach, S. Szymczak, Surrogate minimal depth as an importance measure for variables in random forests, Bioinformatics 2019, 35, 3663-3671.


# Install
```
library(devtools)
install_github("StephanSeifert/SurrogateMinimalDepth")
```

# Example data
The package contains an example data set which consists of a single replicate of the simulation study 1 in our manuscript. Please refer to the paper and the documentation of the SMD_example_data for further details on the simulation scenario.

# Usage
First the package and the example data are loaded:
```
library(SurrogateMinimalDepth)
data("SMD_example_data")
dim(SMD_example_data)
[1] 100 201

head(SMD_example_data[, 1:5])
           y          X1         X2         X3           X4
1  1.8222421 -0.02768266 -1.1019154  2.2659401  0.008021516
2 -1.0401813  0.73258486 -0.4107975  0.7587792 -0.718752746
3  2.7139607 -0.05399936  1.1851261  0.9743160 -2.563176970
4 -0.7081372 -0.84838121 -0.8975802  0.5247899  1.180683275
5 -1.0264429 -0.42219003  0.5439467 -0.1626504  0.682333020
6  3.1871209  0.91722598  0.1974106  0.9571554  0.351634641

```
The data set has 100 observations in the rows and the columns contain the continuous outcome variable y and 200 continuous predictor variables in the columns.

## Minimal depth

First, we perform variable selecion based on minimal depth using 1000 trees in the random forest. To make the analysis reproducible we set the seed first.
```
set.seed(42)
res.md = var.select.md(x = SMD_example_data[,2:ncol(SMD_example_data)], y = SMD_example_data[,1], ntree=1000)
res.md$var
[1] "X2"     "X3"     "X4"     "X5"     "X6"     "cp1_8"  "cp2_6"  "cp2_7"  "cp3_4"  "cp3_6"  "cp8_10" "cgn_68" "cgn_72" "cgn_81"
```

The selected variables are stored in res.md$var. In this analysis the relevant basic variables X2 to X6, as well as the relevant variables cp2_6, cp2_7, cp3_4, and cp3_6, and the non-relevant variables cp8_10, cgn_68, cgn_72, and cgn_81 are selected.

The MD values for each predictor variable and the threshold to select variables can be extracted as follows:
```
md = res.md$info$depth
head(md)
   X1    X2    X3    X4    X5    X6 
9.823 7.848 6.164 6.662 6.442 6.390 

res$info$threshold
[1] 9.23097
```
We can see that variables X2, …, X6 have MD values smaller than the threshold in contrast to X1.

## Surrogate Minimal depth
Now we would like to analyze the example data with surrogate minimal depth which works similarly. However, we need to specify an additional parameter s, i.e. the number of surrogate variables that should be considered. In this analysis we use s = 10. Based on our simulation studies we recommend to set this parameter to approximately 1% of the predictor variables in larger datasets. 

Variable selection with var.select.smd is conducted:
```
set.seed(42)
res.smd = var.select.smd(x = SMD_example_data[,2:ncol(SMD_example_data)], y = SMD_example_data[,1], s = 10, ntree = 1000)
res.smd$var
 [1] "X1"     "X2"     "X3"     "X4"     "X5"     "X6"     "cp1_1"  "cp1_2"  "cp1_3"  "cp1_4"  "cp1_5"  "cp1_6"  "cp1_7"  "cp1_8"  "cp1_9" 
[16] "cp1_10" "cp2_4"  "cp2_6" 
```


The selected variables are stored in res.smd$var. In this analysis the relevant basic variables X1 to X6, as well as the relevant variables cp1_1 to cp1_10, cp2_4, and cp2_6 are selected. Compared to MD more of the relevant variables and none of the non-relevant variables are selected.


The SMD values for each predictor variable and the threshold to select variables can be extracted as follows:

```
smd = res.smd$info$depth
head(smd)
  X1    X2    X3    X4    X5    X6 
2.344 2.287 2.095 2.576 2.509 2.276 

res.smd$info$threshold
[1] 2.690082
```

We can see that variables X1, …, X6 have SMD values smaller than the threshold.

## Variable relations (based on the mean adjusted agreement of surrogate variables)
Now we want to investigate the relations of variables. We would like to identify which of the first 100 predictor variables are related to X1 and X7. We simulated 10 correlated predictor variables for each of these two basic variables.
One possibility to investigate variable relations is to use the results from var.select.smd. Hence, first SMD is conducted like in the previous section:

```
res.smd = var.select.smd(x = SMD_example_data[,2:ncol(SMD_example_data)], y = SMD_example_data[,1], s = 10, ntree = 1000)
```
Subsequently, variable relations are analyzed with var.relations. The parameter t can be adapted to either focus on strongly related variables only (high numbers) or to include also moderately related variables (low numbers):

```
candidates = colnames(SMD_example_data )[2:101]
rel = var.relations(forest = res.smd$forest, variables = c("X1","X7"), candidates = candidates, t = 5)
rel$var
$X1
 [1] "cp1_1"  "cp1_2"  "cp1_3"  "cp1_4"  "cp1_5"  "cp1_6"  "cp1_7"  "cp1_8"  "cp1_9"  "cp1_10"

$X7
 [1] "cp7_1"  "cp7_2"  "cp7_3"  "cp7_4"  "cp7_5"  "cp7_6"  "cp7_7"  "cp7_8"  "cp7_9"  "cp7_10"
```

All of the variables that are correlated to X1 are correctly identified as related to X1 and all of the variables that are correlated to X7 are correcly identified as related to X7. 
