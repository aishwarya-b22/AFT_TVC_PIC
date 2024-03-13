# Summarise results and provides plots based on the real dataset:

simSummaryAFT <- function(numPoints, maxTime, estPlots = TRUE, predSurvPlot = TRUE, repNum){
  
  source("fnstore_AFT_TVC_PIC.R")
  source("hessianMat_AFT_TVC_PIC.R")
  library(latex2exp)
  
  if(repNum > repeats){print("Choose a replication number within the simulation study.")}
  
  if(missing(numPoints)){numPoints = 200}
  if(missing(maxTime)){maxTime = 1}
  if(missing(estPlots)){estPlots = TRUE}
  if(missing(predSurvPlot)){predSurvPlot = TRUE}
  if(predSurvPlot == TRUE & missing(repNum)){repNum = 1}
  
  ##############################################
  ### Generating results for beta and gamma: ###
  ##############################################
  
  # Number of replications in simulation study:
  repeats = length(repID)
  
  # Create empty matrix to store results:
  results = matrix(rep(0,15), nrow = 5, ncol = 3)
  colnames(results) = c("Beta_1", "Beta_2", "Gamma")
  rownames(results) = c("Bias", "Monte Carlo SD",
                        "Average Asymptotic SD",
                        "Coverage Probability (MCSD)", "Coverage Probability (AASD)")
  
  numPar = dim(bStore)[2] + dim(gStore)[2] + dim(tStore)[2]
  bgDim = dim(bStore)[2] + dim(gStore)[2]
  bDim = dim(bStore)[2] 
  
  AASDEigenmat = {}
  
  for(i in 1:repeats){
    AASDEigenmat = rbind(AASDEigenmat, sqrt(diag(asyEigenVarStore[[i]])))
  }
  
  allFiniteAASDEigen = which(is.finite(rowSums(AASDEigenmat)))
  percentNaNAASDEigen =  100 - (length(allFiniteAASDEigen)/repeats)*100
  allFinite = allFiniteAASDEigen
  varStore = asyEigenVarStore
  
  # cat(percentNaNAASDEigen, "percent (by reps) of the asymptotic variances has NaN values"\n)
  
  #############################
  # Calculating bias results: #
  #############################
  
  results[1,] = c(apply(bStore[allFinite,], 2, mean), mean(gStore[allFinite])) - c(beta_true, gamma_true)
  
  #############################
  # Calculating MCSD results: #
  #############################
  
  results[2,] = c(sd(bStore[allFinite,1]), sd(bStore[allFinite,2]), sd(gStore[allFinite]))
  
  ####################################################
  # Calculating various asymptotic variance results: #
  ####################################################
  
  results[3,] = apply(AASDEigenmat[allFiniteAASDEigen,1:3], 2, mean)
  
  #######################################
  # Calculating coverage probabilities: #
  #######################################
  b1cvg = 0
  b2cvg = 0
  g1cvg = 0
  
  b1AAcvg = 0
  b2AAcvg = 0
  g1AAcvg = 0
  
  for(i in allFinite){
    
    # Beta_1 Confidence Intervals using MCSD:
    beta1_lower = bStore[i,1] - qnorm(0.975)*results[2,1]
    beta1_upper = bStore[i,1] + qnorm(0.975)*results[2,1]
    
    if(beta_true[1] >= beta1_lower & beta_true[1] <= beta1_upper){b1cvg = b1cvg + 1}
    
    # Beta_2 Confidence Intervals using MCSD:
    beta2_lower = bStore[i,2] - qnorm(0.975)*results[2,2]
    beta2_upper = bStore[i,2] + qnorm(0.975)*results[2,2]
    
    if(beta_true[2] >= beta2_lower & beta_true[2] <= beta2_upper){b2cvg = b2cvg + 1}
    
    # Gamma Confidence Intervals using MCSD:
    gamma_lower = gStore[i] - qnorm(0.975)*results[2,3]
    gamma_upper = gStore[i] + qnorm(0.975)*results[2,3]
    
    if(gamma_true >= gamma_lower & gamma_true <= gamma_upper){g1cvg = g1cvg + 1}
    
    # Beta_1 Confidence Intervals using AASD:
    beta1_lower_AA = bStore[i,1] - qnorm(0.975)*results[3,1]
    beta1_upper_AA = bStore[i,1] + qnorm(0.975)*results[3,1]
    
    if(beta_true[1] >= beta1_lower_AA & beta_true[1] <= beta1_upper_AA){b1AAcvg = b1AAcvg + 1}
    
    # Beta_2 Confidence Intervals using AASD:
    beta2_lower_AA = bStore[i,2] - qnorm(0.975)*results[3,2]
    beta2_upper_AA = bStore[i,2] + qnorm(0.975)*results[3,2]
    
    if(beta_true[2] >= beta2_lower_AA & beta_true[2] <= beta2_upper_AA){b2AAcvg = b2AAcvg + 1}
    
    # Gamma Confidence Intervals using AASD:
    gamma_lower_AA = gStore[i] - qnorm(0.975)*results[3,3]
    gamma_upper_AA = gStore[i] + qnorm(0.975)*results[3,3]
    
    if(gamma_true >= gamma_lower_AA & gamma_true <= gamma_upper_AA){g1AAcvg = g1AAcvg + 1}
  }
  
  results[4,] = c(b1cvg, b2cvg, g1cvg)/length(allFinite)
  
  results[5,] = c(b1AAcvg, b2AAcvg, g1AAcvg)/length(allFinite)
  
  ###########################
  ### Summary of results: ###
  ###########################
  
  print(results)
  
  #################################################################
  ### Generating results for the hazard and survival functions: ###
  #################################################################
  
  if(estPlots == TRUE){
    
    outKappa = numPoints
    maxVal = maxTime
    kappaGrid = seq(0, maxVal, length.out = outKappa)
    
    hazMat = {}
    survMat = {}
    
    for(i in 1:repeats){
      gaussObj = gauss_basis_fun(tStore[i,], kappaGrid, knotStore[i,], sdStore[i,])
      
      hazard = gaussObj$haz_gen
      hazMat = rbind(hazMat, hazard)
      rownames(hazMat) = NULL
      
      cumHazard = gaussObj$c_haz_gen
      survVec = exp(-cumHazard)
      survMat = rbind(survMat, survVec)
      rownames(survMat) = NULL
    }
    
    hazMatSort = apply(hazMat, 2, sort)
    
    avgHaz = apply(hazMatSort, 2, mean)
    avgSurv = apply(survMat, 2, mean)
    
    #################################################################
    ### Plot hazard/survival functions with confidence intervals: ###
    #################################################################
    
    if(current_dist == "weibull"){
      trueHazVals = alpha*psi*kappaGrid^{alpha-1}
      trueSurvVals = exp(-psi*kappaGrid^alpha)
    }else if(current_dist == "log-logistic"){
      trueHazVals = (alpha*(kappaGrid^{alpha-1})*psi)/(1+(psi*kappaGrid^{alpha}))
      trueSurvVals = 1/(1 + psi*kappaGrid^alpha)
    }
    
    # Plotting average confidence interval bounds
    lowerHazInt = {}
    upperHazInt = {}
    
    lowerSurvInt = {}
    upperSurvInt = {}
    
    for(i in allFinite){
      dataHazVar = {}
      dataSurvVar = {}
      
      hazVec = hazMat[i,]
      survVector = survMat[i,]
      
      varTheta = varStore[[i]][(bgDim+1):numPar, (bgDim+1):numPar]
      
      for(j in 1:outKappa){
        gaussObj = gauss_basis_fun(tStore[i,], kappaGrid[j], knotStore[i,], sdStore[i,])
        
        gsb = gaussObj$gsb
        varHaz = t(gsb)%*%varTheta%*%gsb
        dataHazVar = c(dataHazVar, varHaz)
        
        cumgsb = gaussObj$cumgsb
        mu = gaussObj$c_haz_gen
        sigsq = t(cumgsb)%*%varTheta%*%cumgsb
        varSurv = exp(-2*mu + sigsq)*(exp(sigsq) - 1)
        dataSurvVar = c(dataSurvVar, varSurv)
      }
      
      lowerHazInt = rbind(lowerHazInt, hazVec - qnorm(0.975)*sqrt(dataHazVar))
      upperHazInt = rbind(upperHazInt, hazVec + qnorm(0.975)*sqrt(dataHazVar))
      
      lowerSurvInt = rbind(lowerSurvInt, survVector - qnorm(0.975)*sqrt(dataSurvVar))
      upperSurvInt = rbind(upperSurvInt, survVector + qnorm(0.975)*sqrt(dataSurvVar))
    }
    
    avgLowerHazCIBound = apply(lowerHazInt, 2, mean)
    avgUpperHazCIBound = apply(upperHazInt, 2, mean)
    
    avgLowerSurvCIBound = apply(lowerSurvInt, 2, mean)
    avgUpperSurvCIBound = apply(upperSurvInt, 2, mean)
    
    #max_y1 = max(avgUpperHazCIBound) + 0.5
    max_y1 = 7
    
    plot(kappaGrid, avgHaz, ylim = c(0, max_y1), type = "l", lwd = 1.5, 
         xlab = expression(kappa), ylab = (TeX(r'($\lambda_0(\kappa)$)')), col = "white")
    grid(lty = "solid", col = "lightgray")
    lines(kappaGrid, avgLowerHazCIBound, col = "#53868B", type = "l", lty = 3, lwd = 2)
    lines(kappaGrid, avgUpperHazCIBound, col = "#53868B", type = "l", lty = 3, lwd = 2)
    polygon(c(kappaGrid, rev(kappaGrid)), c(avgLowerHazCIBound, rev(avgUpperHazCIBound)), col = rgb(0.56, 0.9, 0.93, alpha = 0.85), border = NA)
    lines(kappaGrid, avgHaz, ylim = c(0, max_y1), type = "l", lty =2, lwd = 2)
    lines(kappaGrid, trueHazVals, ylim = c(0, max_y1), type = "l", col = "red", lwd = 2)
    legend(x = "topleft",          # Position
           legend = c("True", "Mean Estimate", "95% AASD CI"),  # Legend texts
           lty = c(1, 2, 3),           # Line types
           col = c("red", "black", "#53868B"),           # Line colors
           lwd = 2,
           cex = 0.9)
    
    plot(kappaGrid, avgSurv, ylim = c(0, 1), xlim = c(0, 1.5), type = "l", lwd = 1.5, 
         xlab = expression(kappa), ylab = (TeX(r'($S_0(\kappa)$)')), col = "white")
    grid(lty = "solid", col = "lightgray")
    lines(kappaGrid, avgLowerSurvCIBound, col = "#53868B", type = "l", lty = 3, lwd = 2)
    lines(kappaGrid, avgUpperSurvCIBound, col = "#53868B", type = "l", lty = 3, lwd = 2)
    polygon(c(kappaGrid, rev(kappaGrid)), c(avgLowerSurvCIBound, rev(avgUpperSurvCIBound)), col = rgb(0.56, 0.9, 0.93, alpha = 0.85), border = NA)
    lines(kappaGrid, avgSurv, ylim = c(0,1), type = "l", lty =2, lwd = 2)
    lines(kappaGrid, trueSurvVals, ylim = c(0,1), type = "l", col = "red", lwd = 2)
    legend(x = "bottomleft",          # Position
           legend = c("True", "Mean Estimate", "95% AASD CI"),  # Legend texts
           lty = c(1, 2, 3),           # Line types
           col = c("red", "black", "#53868B"),           # Line colors
           lwd = 2,
           cex = 0.9) 
  }
  
  if(predSurvPlot == TRUE){
    
    source("dataGen_AFT_TVC_PIC.R")
    set.seed(repNum + 100)
    dataPred = dataGenAFT(n, beta_true, gamma_true, tau_min, tau_max, current_dist, alpha, psi, pi_E = event_prop, alpha_L, alpha_R)
    
    betaVal = bStore[i,]
    gammaVal = gStore[i]
    thetaVal = tStore[i,]
    asyVar = asyEigenVarStore[[i]]
    basisKnots = knotStore[i,]
    basisSD = sdStore[i,]
    activeTheta = as.numeric(abs(thetaVal) < 0.01 & (allGrad[i,])[(bgDim + 1) :numPar] < 0)
    
    Xmat = dataPred$Xmat
    delE = dataPred$del[[1]]
    delL = dataPred$del[[2]]
    delR = dataPred$del[[3]]
    delI = dataPred$del[[4]]
    censor_type = delE*1 + delL*3 + delI*4 + delR*2
    yR = dataPred$yR
    yL = dataPred$yL
    t_ni = (yR^delE)*(yR^delL)*(yL^delR)*(yR^delI)
    max_time = max(t_ni) 
    
    data = list(Xmat = cbind(1:nrow(Xmat), Xmat), censor_type = censor_type, max_time = max_time)
    
    std = sqrt(diag(asyVar))
    
    estimates_AB_simulation = list(beta_hat = betaVal, gamma_hat = gammaVal, theta_hat = thetaVal, asym_std_beta = std[1:bDim], asym_std_gamma = std[bgDim], asym_std_theta = std[c((bgDim +1):numPar)[which(activeTheta != 1)]], knots_location = basisKnots, knots_sigma = basisSD, cov_mat_theta = asyVar[(bgDim +1):numPar, (bgDim +1):numPar], index_active_theta = activeTheta)
    rm(std)
    
    source("plAFT_GIC_TVC_Gaussian_1_autosmooth_AB_simulation_plot_design.R")
    
    plAFT_GIC_TVC_Gaussian_1_autosmooth_AB_simulation_plot_design(data = data, estimation_result = estimates_AB_simulation, knots_option = "percentile", plot_sig_lv = 0.05)
    
  }
  
}