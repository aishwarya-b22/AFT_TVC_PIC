---
Title: "Simulation and Inference"
Author: "Aishwarya Bhaskaran"
Date: "March 2024"
Output:
  html_document: default
  pdf_document: default
---

## Data Generation

In this section, we briefly provide details on the function used to simulate data from an accelerated failure time model with time-varying covariates and partly interval censoring.
The following R file contains the code to generate such data:

```r
source("dataGen_AFT_TVC_PIC.R")
```

### Input
In order to generate the required data, the following inputs are necessary:

| Input parameter  | Description |
| ------------- | ------------- |
| n | Sample size  |
| beta   | True values of beta  |
| gamma  | True values of gamma  |
| tau_min   | To generate treatment times according to Unif(tau_min, tau_max). The default value is 0. |
| tau_max  | To generate treatment times according to Unif(tau_min, tau_max). The default value is 1. |
| dist  | Distribution of hazard; choice between "weibull" and "log-logistic" |
| alpha  | Shape parameter for specifying the distribution  |
| psi  | Scale parameter for specifying the distribution  |
| pi_E  | Average proportion of events in simulated sample  |
| alpha_L   | Decreasing this value lowers proportion of left-censoring & increases proportion of interval censoring. The default value is 0.5.   |
| alpha_R  | Increasing this value lowers proportion of right-censoring & increases proportion of interval censoring. The default value is 1.5.  |

### Simulating data from an AFT model with time-varying covariates and partly-interval censoring
Here is an example that uses the R function to simulate the required data.
```r
data = dataGenAFT(n = 100, beta = c(1, -1), gamma = -0.5, tau_min = 0, tau_max = 2, dist = "weibull", alpha = 3, psi = 1, pi_E = 0.3, alpha_L = 0.9, alpha_R = 1.1)
```
### Output
The following list of outputs are useful and are necessary (except data$exact) for the optimisation process in the next section.

| Output parameter  | Description |
| ------------- | ------------- |
| exact | Failure times for all individuals ($T_i, 1 \leq i \leq n$) |
| yL  | Left time points ($y_i^L, 1 \leq i \leq n$)|
| yR  | Right time points ($y_i^R, 1 \leq i \leq n$)|
| Xmat  | Matrix with time-fixed covariates  |
| Zmat   | Matrix with time-varying covariates (Long-format)  |
| zfin | $z_i(t_{ni})$; values of time-varying covariates at $t_{ni}$)
| del  | List of all censoring indicators (delE ($\delta_i$), delL ($\delta_i^L$), delR ($\delta_i^R$), delI ($\delta_i^I$)) |
| censor_prop | Proportions of each type of censoring |

Here we can preview the censoring proportions in the simulated dataset.
```r
attach(data)
print(censor_prop)
```
```
      Event times     Left-censored    Right-censored Interval-censored 
             0.27              0.45              0.12              0.16
```

We can preview the short-format and long-format data below.

Short-format data:
```r
short = head(cbind(data$yL, data$yR, Xmat, zfin, del[[1]], del[[2]], del[[3]], del[[4]]))
colnames(short) = c("yL", "yR", "x1", "x2", "z_i(t_ni)", "delE", "delL", "delR", "delI")
head(short)
```
```
             yL         yR x1       x2 z_i(t_ni) delE delL delR delI
[1,] 0.12537003 0.79929490  1 1.374922         1    0    0    0    1
[2,] 0.14411819 0.14411819  0 1.240363         0    1    0    0    0
[3,] 0.00000000 0.79772737  1 2.863321         1    0    1    0    0
[4,] 0.09430258 0.09430258  0 2.341156         0    1    0    0    0
[5,] 0.18171990 0.18171990  0 1.537551         0    1    0    0    0
[6,] 0.00000000 0.37218399  0 2.547797         0    0    1    0    0
```

Long-format data:
```r
head(Zmat)
```
```
     Individual      Start        End z_i(t)
[1,]          1 0.00000000 0.08972662      0
[2,]          1 0.08972662 0.79929490      1
[3,]          2 0.00000000 0.14411819      0
[4,]          3 0.00000000 0.43621940      0
[5,]          3 0.43621940 0.79772737      1
[6,]          4 0.00000000 0.09430258      0
```

## Optimisation
In this section, we briefly provide details on the optimisation function used to obtain estimates of the regression and spline coefficients and as well as the baseline hazard, via a maximum penalised likelihood approach. The R file below contains the code to optimise the parameters. It can be used on either simulated data or real-world data.

```r
source("optim_AFT_TVC_PIC.R")
```
### Input
We provide descriptions of any additional inputs required for the optimisation function below:

| Input parameter  | Description |
| ------------- | ------------- |
| m | Number of basis functions to be used  |
| knots_option | Choice of either "quantile" or "equal-space" to select the locations of the basis functions |
| knotSTOP | Maximum number of iterations afterwhich the basis functions are fixed. The default value is 100 |
| sd_option | Choice of either 1 (Standard deviation is half of widths between knots) or 2 (Standard deviation chosen ensures that 2/3 of the data points are covered by each basis function). The default value is 1 |
| quantVec | Probabilities of the quantiles associated with the first and last basis functions. The default value is c(0.05, 0.95) |
| h_init | Starting value for the smoothing parameter. The default value is 1 |
| smooth_stop | Set as FALSE if automatic smoothing selection is required. The default value is FALSE |
| tol1 | Tolerance level for change in estimates. The default value is 1e-4 |
| tol2 | Tolerance level for gradients of the regression cofficients. The default value is 1e-4 |
| maxIterPerLoop | Maximum number of iterations per loop. The default value is 1000 |
| diffDF | Tolerance for change in degrees of freedom (to establish convergence of the smoothing parameter). The default value is 1 |
| stableNum | Number of consecutive outer loops that diffDF has to be satisfied for. The default value is 1 |
| outerMax | Maximum number of outer loops. The default value is 10 |

### Optimising a semiparametric AFT model with time-varying covariates and partly-interval censoring
Here is a code snippet that uses the R function optimAFT( ) with sample inputs to optimise the semiparametric AFT model:
```r
optObj = optimAFT(n = 100, m = 5, Xmat = data$Xmat, Zmat = data$Zmat, zfin = data$zfin, yL = data$yL, yR = data$yR, delE = data$del[[1]], delL = data$del[[2], delR = data$del[[3]], delI = data$del[[4]], knots_option = "quantile", knotSTOP = 100, sd_option = 2 , quantVec = c(0.05, 0.95), h_init = 1, smooth_stop = FALSE, tol1 = 1e-6, tol2 = 1e-6, maxIterPerLoop = 5000, diffDF = 0.5, stableNum = 3, outerMax = 10)
```
### Output
The following list is a list of useful outputs from the optimisation function above.

| Output parameter  | Description |
| ------------- | ------------- |
| bVal | Estimates of beta (regression coefficients for time-fixed covariates) | 
| gVal | Estimates of gamma (regression coefficients for time-varying covariates) | 
| tVal | Estimates of theta (spline coefficients) | 
| grad | Gradients of all model parameters at convergence | 
| h | Final value of the smoothing parameter |
| iter | Total number of iterations taken for the algorithm to converge |
| count | Total number of outer loops required for the algorithm to converge | 
| rep_knots | Final locations used for the basis functions| 
| rep_sd | Final standard deviations used for the basis functions |
| haz | Estimated baseline hazard at the point of $\kappa$ |
| asyEigenVar | Asymptotic variance-covariance matrix of the model parameters at convergence |
| cvg | 1 if the algorithm has converged and 0 otherwise |  
| logLikVec | Vector of the penalised log-likelihood values recorded at the end of every iteration |

We will present examples of the output (for several replicates) and discuss them in greater details in the next few sections.

## Simulation
In our simulation, we perform a pre-determined number of replications that replicate both the data generation and optimisation process to obtain several sets of estimates. The following R file runs simulation studies bases on the input specified.

```r
source("simStudy_AFT_TVC_PIC.R")
```

### Input
In addition to the input parameters mentioned earlier, one needs to specify the number of repeats required for the simulation study:

| Input parameter  | Description |
| ------------- | ------------- |
| repeats | Number of replicates desired for the simulation study  |

### Simulation Study
Here is a code snippet that uses the R function simRunAFT( ) with sample inputs to run a simulation study with 200 replicates based on a Weibull hazard with a sample size of $n = 100$ each. The average event proportion was set as 0.3 and $m$ (number of basis functions) was set to be 5.:
```r
simRes = simRunAFT(repeats = 200, n = 1000, m = 6, beta_true  = c(1, -1), gamma_true = -0.5, event_prop = 0.7, current_dist = "weibull", alpha = 3, psi = 1,
          tau_min = 0, tau_max = 2, alpha_L = 0.9, alpha_R = 1.1, knotsOpt = "quantile", knotPlace = 2, quantVec = c(0.29, 0.99), sdOpt = 2, 
          knotMaxIter  = 100, h_init = 0.001, smooth_stop = FALSE, maxIterPerLoop = 2000, tol1 = 1E-4, tol2 = 1E-4, diffDF = 0.7, stableNum = 1, outerMax = 6)
```
We summarise the results and perform inference in the next section.

## Summarising Simulation Results
In this section, we summarise the results, provide coverage probabilities and produce plots based on the simulation results. We first load the data file containing the results from the simulation study conducted in the previous section. With this, we can also replicate the plots in Figure 3 of the manuscript.

```r
load("AFT_TVC_PIC_weibull_E0.7_n1000_m10.RData")
```
The following R file runs the summary functions for the results obtained from the simulation study.

```r
source("simSummary_AFT_TVC_PIC.R")
```

### Input
In order to use the summary function for the simulation results, simulSummaryAFT( ), one needs to specify the following inputs:

| Input parameter  | Description |
| ------------- | ------------- |
| numPoints | Number of points in the grid of accelerated failure times to generate the estimated baseline hazard and survival plots |
| maxTime | Specifies the maximum value for the range for the grid of accelerated failure times for the estimated baseline hazard and survival plots |
| estPlots | If TRUE, the estimated baseline hazard and survival plots will be generated |
| predSurvPlot | If TRUE, the predictive survival plots will be generated |
| repNum | Replication number from the simulation study for which the predictive survival plots should be generated for. The default value is 1. |

### Results
The following code generates the values of the bias, Monte Carlo standard deviation (MCSD), average asymptotic standard deviation (AASD) and the coverage probabilities (CP) calculated from generating 95\% confidence intervals using both the MCSDs and AASDs.

```r
simSummaryAFT(numPoints = 200, maxTime = 1.5, estPlots = TRUE, predSurvPlot = TRUE, repNum = 1)
```

```
                                Beta_1       Beta_2        Gamma
Bias                        0.00102714 -0.008103535 -0.003022216
Monte Carlo SD              0.02801281  0.048015376  0.038969996
Average Asymptotic SD       0.02696118  0.041693834  0.038637033
Coverage Probability (MCSD) 0.95500000  0.950000000  0.960000000
Coverage Probability (AASD) 0.94000000  0.915000000  0.960000000
```
The following estimated baseline hazard and survival plots are also generated:

[github_weibull_E0.7_n1000_hazard.pdf](https://github.com/aishwarya-b22/AFT_TVC_PIC/files/14584350/github_weibull_E0.7_n1000_hazard.pdf)

## Summarising Application Results

```r
source("realDataSummary_AFT_TVC_PIC.R")
``` 

```r
data = list(Xmat = read.csv("WBRT_Xmat.csv")[, -1], tmat = read.csv("WBRT_tmat_PIC.csv")[, -1])


realDataSummaryAFT(numPoints = 200, maxTime = quantile(postOpt$kappa_vec, 0.75), estPlots = TRUE, predSurvPlot = TRUE)
```
