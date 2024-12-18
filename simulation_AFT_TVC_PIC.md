---
Title: "Simulation and Inference"
Author: "Aishwarya Bhaskaran"
Date: "December 2024"
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
| tau_min  | To generate treatment times according to Unif(tau_min, tau_max). The default value is 0. |
| tau_max  | To generate treatment times according to Unif(tau_min, tau_max). The default value is 1. |
| dist  | Distribution of hazard; choice between "weibull" and "log-logistic" |
| alpha  | Shape parameter for specifying the distribution  |
| psi  | Scale parameter for specifying the distribution  |
| pi_E  | Average proportion of events in simulated sample  |
| alpha_L  | Decreasing this value lowers proportion of left-censoring & increases proportion of interval censoring. The default value is 0.5. |
| alpha_R  | Increasing this value lowers proportion of right-censoring & increases proportion of interval censoring. The default value is 1.5. |

### Simulating data from an AFT model with time-varying covariates and partly-interval censoring
Here is an example that uses the R function to simulate the required data.
```r
set.seed(123)
data = dataGenAFT(n = 100, beta = c(1, -1), gamma = -0.1, tau_min = 0, tau_max = 2, dist = "weibull", alpha = 3, psi = 1, pi_E = 0.3, alpha_L = 0.9, alpha_R = 1.1)
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
             0.33              0.41              0.14              0.12
```

We can preview the short-format and long-format data below.

Short-format data:
```r
short = head(cbind(data$yL, data$yR, Xmat, zfin, del[[1]], del[[2]], del[[3]], del[[4]]))
colnames(short) = c("yL", "yR", "x1", "x2", "z_i(t_ni)", "delE", "delL", "delR", "delI")
head(short)
```
```
             yL        yR x1        x2 z_i(t_ni) delE delL delR delI
[1,] 0.00000000 0.3182455  0 1.7999669         0    0    1    0    0
[2,] 0.30744422 0.3074442  1 0.9984706         1    1    0    0    0
[3,] 0.00000000 0.2583901  0 1.4658391         0    0    1    0    0
[4,] 0.07197562 0.4102962  1 2.8634215         0    0    0    0    1
[5,] 0.52344045       Inf  1 1.4487072         0    0    0    1    0
[6,] 0.00000000 0.1602124  0 2.6710507         0    0    1    0    0
```

Long-format data:
```r
head(Zmat)
```
```
     Individual      Start        End z_i(t)
[1,]          1 0.00000000 0.31824547      0
[2,]          2 0.00000000 0.01885981      0
[3,]          2 0.01885981 0.30744422      1
[4,]          3 0.00000000 0.25839012      0
[5,]          4 0.00000000 0.41029622      0
[6,]          5 0.00000000 0.52344045      0
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
simRes = simRunAFT(repeats = 200, n = 100, m = 5, beta_true  = c(1, -1), gamma_true = -0.1, event_prop = 0.3, current_dist = "weibull", alpha = 3, psi = 1,
          tau_min = 0, tau_max = 2, alpha_L = 0.9, alpha_R = 1.1, knotsOpt = "quantile", knotPlace = 1, quantVec = c(0.1, 0.9), sdOpt = 2, 
          knotMaxIter  = 100, h_init = 0.001, smooth_stop = FALSE, maxIterPerLoop = 2000, tol1 = 1E-6, tol2 = 1E-6, diffDF = 0.5, stableNum = 3, outerMax = 10)
```
We summarise the results and perform inference in the next section.

## Summarising Simulation Results
In this section, we summarise the results, provide coverage probabilities and produce plots based on the simulation results. We first load the data file containing the results from the simulation study conducted in the previous section. With this, we can also replicate some of the plots in the manuscript.

```r
load("AFT_TVC_PIC_weibull_E0.3_n1000_m10.RData")
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
                                Beta_1      Beta_2      Gamma
Bias                        -0.0036853 -0.05452317 -0.1100646
Monte Carlo SD               0.1087588  0.07918904  0.2171554
Average Asymptotic SD        0.1174895  0.06124770  0.1932592
Coverage Probability (MCSD)  0.9500000  0.89000000  0.9250000
Coverage Probability (AASD)  0.9650000  0.80500000  0.8650000
```
The following estimated baseline hazard and survival plots are generated:

[weibull_haz_E0.3_n100_m5.pdf](https://github.com/user-attachments/files/18174539/weibull_haz_E0.3_n100_m5.pdf)
[weibull_surv_E0.3_n100_m5.pdf](https://github.com/user-attachments/files/18174540/weibull_surv_E0.3_n100_m5.pdf)

## Summarising Application Results
In this section, we summarise the results and produce dynamic prediction plots based on the results of the WBRTMel trial dataset. We first load the data file containing the results from optimising the real dataset. With this, we can also replicate some of the plots in the manuscript.

```r
load("realData_systemictherapy.RData")
```
The following R file runs the summary functions for the results obtained from the analysis of the dataset used in the application.

```r
source("realDataSummary_AFT_TVC_PIC.R")
```

### Results
The following code generates the estimated coefficients and as well as the predictive survival plots for the various time-fixed and time-varying covariates.

```r
data = list(Xmat = read.csv("WBRT_Xmat.csv")[, -1], tmat = read.csv("WBRT_tmat_PIC.csv")[, -1])
realDataSummaryAFT(numPoints = 200, maxTime = quantile(postOpt$kappa_vec, 0.75), estPlots = FALSE, predSurvPlot = TRUE)
```
```
     Estimated coefficient pValue                 sigResult
[1,] "0.574236056606409"   "0.000161588412554218" "***"    
[2,] "-0.472578010930059"  "0.000978836551845801" "***"    
[3,] "-0.704817699590601"  "1.75378204919216e-06" "***"    
[4,] "-0.0365240999833249" "0.794999109528828"    "NSF"    
[5,] "-0.0254496441606496" "1.99929005111838e-28" "***"    
[6,] "0.482898201981513"   "0.0492377260867074"   "*"
```
We can also generate predictive survival plots. As an example, here we present a concise version of Figure 5 that displays plots of predicted distant intracranial free survival
(in years) for a time-fixed covariate of interest, treatment (WBRT vs observation).

![github_WBRT_predSurv_treatment](https://github.com/aishwarya-b22/AFT_TVC_PIC/assets/61529713/c2056081-52cf-4e96-9154-9fb647266322)


