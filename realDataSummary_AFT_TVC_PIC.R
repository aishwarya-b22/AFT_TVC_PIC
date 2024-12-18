realDataSummaryAFT <- function(numPoints, maxTime, estPlots = TRUE, predSurvPlot = TRUE){
  
  source("fnstore_AFT_TVC_PIC.R")
  source("hessianMat_AFT_TVC_PIC.R")
  library(latex2exp)
  
  if(missing(numPoints)){numPoints = 200}
  if(missing(maxTime)){maxTime = 1}
  if(missing(estPlots)){estPlots = TRUE}
  if(missing(predSurvPlot)){predSurvPlot = TRUE}
  
  # Estimates of Beta, Gamma and Theta:
  bVal = postOpt$bVal
  gVal = postOpt$gVal
  tVal = postOpt$tVal
  
  bgDim = length(bVal) + length(gVal)
  numPar = bgDim + length(tVal)
  bDim = dim(bVal)[2] 
  
  activeTheta = as.numeric(abs(tVal) < 0.01 & postOpt$grad[(bgDim+1):numPar] < 0)
  basisKnots = postOpt$rep_knots
  basisSD = postOpt$rep_sd
  
  # Asymptotic Standard Deviation:
  asyVar = postOpt$asyEigenVar
  asySD = sqrt(diag(asyVar))
  
  # Confidence Intervals for Beta and Gamma:
  confLower = c(bVal, gVal) - qnorm(0.975)*asySD[1:bgDim]
  confUpper = c(bVal, gVal) + qnorm(0.975)*asySD[1:bgDim]
  
  # Wald's test results:
  waldStats = c(bVal, gVal)/asySD[1:bgDim]
  pValue = (pnorm(abs(waldStats), lower.tail = FALSE))*2
  
  sigResults = rep(0, bgDim)
  sigResults[pValue <= 1] = "NSF"
  sigResults[pValue <= 0.1] = "."
  sigResults[pValue <= 0.05] = "*"
  sigResults[pValue <= 0.01] = "**"
  sigResults[pValue <= 0.001] = "***"
  
  # Results:
  res = cbind(c(bVal, gVal), pValue, sigResults)
  colnames(res) = c("Estimated coefficient", "pValue", "sigResult")
  print(res)
  
  # Estimated baseline hazard and baseline survival:
  outKappa = numPoints
  maxVal = maxTime
  kappaGrid = seq(0, maxVal, length.out = outKappa)
  
  gaussObj = gauss_basis_fun(tVal, kappaGrid, postOpt$rep_knots, postOpt$rep_sd)
  
  hazVec = gaussObj$haz_gen
  cumHazard = gaussObj$c_haz_gen
  survVec = exp(-cumHazard)
  
  # Plotting average confidence interval bounds
  if(estPlots == TRUE){
    lowerHazInt = {}
    upperHazInt = {}
    
    lowerSurvInt = {}
    upperSurvInt = {}
    
    dataHazVar = {}
    dataSurvVar = {}
    
    varTheta = asyVar[(bgDim+1):numPar, (bgDim+1):numPar]
    
    for(j in 1:outKappa){
      gaussObj = gauss_basis_fun(tVal, kappaGrid[j], postOpt$rep_knots, postOpt$rep_sd)
      
      gsb = gaussObj$gsb
      varHaz = t(gsb)%*%varTheta%*%gsb
      dataHazVar = c(dataHazVar, varHaz)
      
      cumgsb = gaussObj$cumgsb
      mu = gaussObj$c_haz_gen
      sigsq = t(cumgsb)%*%varTheta%*%cumgsb
      varSurv = exp(-2*mu + sigsq)*(exp(sigsq) - 1)
      dataSurvVar = c(dataSurvVar, varSurv)
    }
    
    lowerHazInt = hazVec - qnorm(0.975)*sqrt(dataHazVar)
    upperHazInt = hazVec + qnorm(0.975)*sqrt(dataHazVar)
    
    lowerSurvInt = survVec - qnorm(0.975)*sqrt(dataSurvVar)
    upperSurvInt = survVec + qnorm(0.975)*sqrt(dataSurvVar)
    
    plot(kappaGrid, hazVec, type = "l", ylim = c(0, max(upperHazInt)))
    lines(kappaGrid, lowerHazInt, col = "red")
    lines(kappaGrid, upperHazInt, col = "red")
    polygon(c(kappaGrid, rev(kappaGrid)), c(lowerHazInt, rev(upperHazInt)), col = rgb(1,0,0, alpha = 0.15), border = NA)
    
    plot(kappaGrid, survVec, type = "l")
    lines(kappaGrid, lowerSurvInt, col = "blue")
    lines(kappaGrid, upperSurvInt, col = "blue")
    polygon(c(kappaGrid, rev(kappaGrid)), c(lowerSurvInt, rev(upperSurvInt)), col = rgb(0,0,1, alpha = 0.15), border = NA)
  }
  
  # Plotting predictive survival curves:
  if(predSurvPlot == TRUE){
    
    std = sqrt(diag(asyVar))
    estimates = list(beta_hat = bVal, gamma_hat = gVal, theta_hat = tVal, asym_std_beta = std[1:5], asym_std_gamma = std[6], asym_std_theta = std[c(7:11)[which(activeTheta != 1)]], knots_location = basisKnots, knots_sigma = basisSD, cov_mat_theta = asyVar[7:11, 7:11], index_active_theta = activeTheta)
    rm(std)
    
    source("plAFT_GIC_TVC_Gaussian_1_autosmooth_WBRT_plot.R")

    plAFT_GIC_TVC_Gaussian_1_autosmooth_WBRT_plot(data = data, estimation_result = estimates, tau = 0)
    plAFT_GIC_TVC_Gaussian_1_autosmooth_WBRT_plot(data = data, estimation_result = estimates, tau = 0.5)
    plAFT_GIC_TVC_Gaussian_1_autosmooth_WBRT_plot(data = data, estimation_result = estimates, tau = 1)
    plAFT_GIC_TVC_Gaussian_1_autosmooth_WBRT_plot(data = data, estimation_result = estimates, tau = 1.5)
    
  }
}


