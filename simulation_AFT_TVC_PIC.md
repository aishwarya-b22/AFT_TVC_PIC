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
| tau_min   | To generate treatment times according to Unif(tau_min, tau_max)  |
| tau_max  | To generate treatment times according to Unif(tau_min, tau_max)  |
| dist  | Distribution of hazard; choice between "weibull" and "log-logistic"  |
| alpha  | Shape parameter for specifying the distribution  |
| psi  | Scale parameter for specifying the distribution  |
| pi_E  | Average proportion of events in simulated sample  |
| alpha_L   | Decreasing this value lowers proportion of left-censoring & increases proportion of interval censoring  |
| alpha_R  | Increasing this value lowers proportion of right-censoring & increases proportion of interval censoring  |

### Simulating data from an AFT model with time-varying covariates and partly-interval censoring
Here is an example that uses the R function to simulate the required data.
```r
data = datagen_AFT_PIC_TVC(n = 100, beta = c(1, -1), gamma = -0.5, tau_min = 0, tau_max = 2, dist = "weibull", alpha = 3, psi = 1, pi_E = 0.3, alpha_L = 0.9, alpha_R = 1.1)
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

## Simulation

## Summarising Simulation Results

## Summarising Application Results
