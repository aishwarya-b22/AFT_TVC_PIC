# Compile x repeats of data simulation and optimisation:
source("sim_AFT_TVC_PIC.R")
source("optim_AFT_TVC_PIC.R")

simRunAFT <- function(repeats, n, m, beta_true, gamma_true, event_prop,
                      current_dist, alpha, psi, tau_min, tau_max, alpha_L, alpha_R,
                      knotsOpt, quantVec, sdOpt, knotMaxIter, h_init, smooth_stop,
                      maxIterPerLoop, tol1, tol2, diffDF, stableNum, outerMax){
  
  # Values to be stored:
  i = 1
  
  repID = {}
  overall_prop = {}
  treat_prop = {}
  
  allIter = {}
  allCount = {}
  
  all_h = {}
  
  smoothRec = list()
  
  bStore = {}
  gStore = {}
  tStore = {}
  
  knotStore = {}
  sdStore = {}
  
  allHaz = {}
  
  asyVarStore = list()
  asyEigenVarStore = list()
  
  allGrad = {}
  logLikStore = list()
  kappaVecStore = list()
  
  allTime = {}
  
  allCvg = {}
  
  start_time <- Sys.time()
  
  while(i <= repeats){
    
    ########################
    #### Simulate data: ####
    ########################
    
    set.seed(i + 100)
    
    data = dataGenAFT(n = n, beta = beta_true, gamma = gamma_true, 
                               tau_min = tau_min, tau_max = tau_max, dist = current_dist,
                               alpha = alpha, psi = psi, pi_E = event_prop, 
                               alpha_L = alpha_L, alpha_R = alpha_R)
    
    # X matrix and Z matrix (long format):
    Xmat = data$Xmat
    Zmat = data$Zmat
    
    # Values of time-varying covariates -> z_i(t_ni):
    zfin = data$zfin
    
    # Vectors of left and right observed times:
    yL = data$yL
    yR = data$yR
    
    # Extract censoring indicator vectors:
    delE = data$del[[1]]
    delL = data$del[[2]]
    delR = data$del[[3]]
    delI = data$del[[4]]
    
    #################################
    #### Optimization Algorithm: ####
    #################################
    
    cat("Starting replication", i, "now: \n")
    
    opt_start <- Sys.time()
    
    h = h_init
    
    postOpt = optimAFT(n = n, m = m, Xmat = Xmat, Zmat = Zmat, zfin = zfin, 
                              yL = yL, yR = yR, delE = delE, delL = delL, delR= delR,
                              delI = delI, knots_option = knotsOpt,
                              knotSTOP = knotMaxIter, sd_option = sdOpt, quantVec = quantVec, h = h, 
                              smooth_stop = smooth_stop, tol1 = tol1, tol2 = tol2,
                              maxIterPerLoop = maxIterPerLoop, diffDF = diffDF, stableNum = stableNum,
                              outerMax = outerMax)
    
    opt_end <- Sys.time()
    
    rep_time = opt_end - opt_start
    
    ###################################################
    ### Outputs to be saved if algorithm converges: ###
    ###################################################
    
    repID = rbind(repID, i)
    
    overall_prop = rbind(overall_prop, data$censor_prop)
    treat_prop = rbind(treat_prop, sum(data$treat_index)/n)
    
    allIter = rbind(allIter, postOpt$iter)
    allCount = rbind(allCount, postOpt$count)
    
    all_h = rbind(all_h, postOpt$h)
    
    h_df_mat = cbind(postOpt$nu_record, postOpt$h_record)
    colnames(h_df_mat) = c("df", "h")
    rownames(h_df_mat) = NULL
    smoothRec = c(smoothRec, list(h_df_mat))
    
    bStore = rbind(bStore, postOpt$bVal)
    gStore = rbind(gStore, postOpt$gVal)
    tStore = rbind(tStore, postOpt$tVal)
    
    knotStore = rbind(knotStore, postOpt$rep_knots)
    sdStore = rbind(sdStore, postOpt$rep_sd)
    
    allHaz = rbind(allHaz, postOpt$haz)
    
    asyEigenVarStore = c(asyEigenVarStore, list(postOpt$asyEigenVar))
    
    allGrad = rbind(allGrad, postOpt$grad)
    
    logLikStore = c(logLikStore, list(postOpt$logLikVec))
    
    kappaVecStore = c(kappaVecStore, list(postOpt$kappa_vec))
    
    allTime = rbind(allTime, rep_time)
    
    
    if(postOpt$cvg == 0){
      rec_error = c(i , postOpt$msg)
      allCvg = rbind(allCvg, rec_error)
    }
    
    cat("Replication", i, "is complete. \n")
    i = i + 1
  }
  
  row.names(repID) = NULL
  end_time <- Sys.time()
  
  total_time = end_time - start_time
  print(total_time)
  
  ###########################################
  ### Save outputs from simulation study: ###
  ###########################################
  
  save(repeats, n, m, beta_true, gamma_true, event_prop, overall_prop, treat_prop,
       current_dist, alpha, psi, repID,  allIter, allCount,
       tau_min, tau_max, alpha_L, alpha_R, knotsOpt, quantVec, sdOpt,
       bStore, gStore, tStore, all_h, knotStore, sdStore, 
       allHaz, asyEigenVarStore, allGrad, logLikStore, 
       allTime, total_time, allCvg, smoothRec, kappaVecStore,
       file = paste0("AFT_TVC_PIC_",current_dist,"_E",event_prop, "_n",n,"_m",m,".RData", collapse = "_"))
  
}