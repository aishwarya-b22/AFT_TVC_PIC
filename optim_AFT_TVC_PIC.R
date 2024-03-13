# Optimization Algorithm:

optimAFT <- function(n, m, Xmat, Zmat, zfin, yL, yR, delE, delL, delR, delI,
                            knots_option, knotSTOP, sd_option, quantVec, h_init,
                            smooth_stop, tol1, tol2, maxIterPerLoop, diffDF, stableNum, outerMax){
  
  if(missing(knotSTOP)){knotSTOP = 100}
  if(missing(sd_option)){sd_option = 1}
  if(missing(quantVec)){quantVec = c(0.05, 0.95)}
  if(missing(h_init)){h_init = 1}
  if(missing(smooth_stop)){smooth_stop = FALSE}
  if(missing(tol1)){tol1 = 1e-4}
  if(missing(tol2)){tol2 = 1e-4}
  if(missing(maxIterPerLoop)){maxIterPerLoop = 1000}
  if(missing(diffDF)){diffDF = 1}
  if(missing(stableNum)){stableNum = 3}
  if(missing(outerMax)){outerMax = 10}
  
  # Set starting values for the parameters to be estimated:
  bVal = rep(0, dim(Xmat)[2])
  gVal = rep(0.5, 1)
  tVal = rep(1, m)
  
  # Required dimensions of vectors of parameters:
  numPar = length(c(bVal, gVal, tVal))
  bgDim = length(c(bVal, gVal))
  
  rep__knots = {}
  rep_sd = {}
  
  Qmat = NA
  Fmat = NA
  FmatEigen = NA
  asyEigenVar = NA
  
  nu_old = -1
  nu_record = -1
  
  h = h_init
  h_record = h_init
  
  oldVal_flag = FALSE
  
  logLikVec = as.vector(0)
  cvg = {}
  msg = {}
  grad = NA
  
  iter = 0
  count = 1
  logLikCount = 0
  nuStableCount = NA
  
  #####################################
  ### Other necessary calculations: ###
  #####################################
  
  # In order to calculate k_i(y_i), i = 1,...n, based on z_i(t_ni):
  id_Z = Zmat[,1]
  time_int = as.matrix(Zmat[,3] - Zmat[,2])
  Z_t = as.matrix(Zmat[,4])
  
  # In order to calculate k_i(y_i^L), for observations with interval censoring:
  count_id_Z = table(id_Z)
  ext_yL = rep(yL, count_id_Z)
  
  # Modified time interval: Only accounts for time intervals up till the left point in interval censoring:
  # First part: Keep entire interval if yL > t_end
  # Second part: Keep partial interval if t_start < yL < t_end
  mod_time_int = ((ext_yL>=Zmat[,3])*(Zmat[,3]-Zmat[,2]) 
                  + (ext_yL>=Zmat[,2] & ext_yL<Zmat[,3])*(ext_yL-Zmat[,2]))
  
  # Indices corresponding to observations that are interval-censored:
  int_ind = which((1:n)*delI != 0)
  
  ########################### Beginning of outer loop: ###########################
  
  repeat{
    iterLoop = 0
    
    cat("Beginning outer loop", count," \n")
    
    ########################### Beginning of inner loop: ###########################
    
    repeat{
      
      iter = iter + 1
      iterLoop = iterLoop + 1
      
      # Calculate vector of kappa values at this iteration:
      kp_right = exp(-Xmat%*%bVal)*as.matrix(tapply(exp(-Z_t%*%gVal)*time_int, id_Z, sum))
      kp_left_int = exp(-Xmat%*%bVal)*as.matrix(tapply(exp(-Z_t%*%gVal)*mod_time_int, id_Z, sum))
      
      if(iter <= knotSTOP){
        
        # Calculate knots and sigma values for Gaussian basis functions:
        kappa_vec = c(kp_right[-int_ind], ((kp_right[int_ind] + kp_left_int[int_ind])/2))
        
        knots_out = knots_fun(num_basis = m, knots_option, quantVec, sd_option, kappa_vec)
        rep_knots = knots_out$gknots
        rep_sd = knots_out$sd
        
        ############################################
        #### Code for R matrix in penalty term: ####
        ############################################
        # Select number of bins between d_1 and d_2:
        num_bins = 5*n
        
        # Find values of bin edges:
        kp_range = range(kappa_vec)
        bin_edges = matrix(seq(kp_range[1],kp_range[2],length.out=num_bins+1))
        
        # Vector of midpoints of bins:
        bin_mid = (bin_edges[1:num_bins] + bin_edges[2:(num_bins+1)])/2
        
        # Calculate R matrix:
        kpstd_bin_mid = apply(as.matrix(bin_mid), 1, function(x){(x-rep_knots)/rep_sd})
        ddpsi_bin_mid = (1/(rep_sd^3))*((kpstd_bin_mid^2)-1)*dnorm(kpstd_bin_mid)
        
        Rmat = (ddpsi_bin_mid%*%t(ddpsi_bin_mid)) * ((kp_range[2]-kp_range[1])/num_bins)
      }
      
      if(iter == knotSTOP + 1){
        tVal[tVal < 0.01] = 0.2
      }
      
      # Calculated updated penalty matrix (h*R):
      h_Rmat = h*Rmat
      
      ########################## Newton step for beta ###########################
      beta_old = bVal
      beta_new = {}
      
      # Baseline hazard and cumulative baseline hazard (most cases):
      set1 = gauss_basis_fun(tVal, kp_right, rep_knots, rep_sd)
      c_haz = set1$c_haz_gen
      haz = set1$haz_gen
      haz_d = set1$haz_d_gen
      
      # Cumulative baseline hazard (When considering the left point in int. censoring):
      set2 = gauss_basis_fun(tVal, kp_left_int, rep_knots, rep_sd)
      c_haz_left = set2$c_haz_gen
      haz_left = set2$haz_gen
      
      #####################################################
      ### Calculate gradient of beta at this iteration: ###
      #####################################################
      expLall = exp(-c_haz_left)*haz_left*kp_left_int
      expRall = exp(-c_haz)*haz*kp_right
      
      ### Fix numerical issues (by adding jitter to the denominator): ###
      
      exp1den = -expm1(-c_haz)  #exp1den = 1-exp(-c_haz)
      exp1den[exp1den == 0] = 1e-6
      
      exp2den = (expm1(-c_haz_left) - expm1(-c_haz))^delI  #exp2den = (exp(-c_haz_left) - exp(-c_haz))^delI
      exp2den[exp2den == 0] = 1e-6
      
      exp3den = haz
      exp3den[exp3den == 0] = 1e-6
      
      # Compute each part in the gradient of beta separately:
      E_part_gb = delE*((-(haz_d*kp_right)/exp3den) -1 + haz*kp_right)
      R_part_gb = delR*(haz*kp_right)
      L_part_gb = delL*(expRall/exp1den)
      I_part_gb = delI*((expLall - expRall)/exp2den)
      
      grad_beta = colSums(as.vector(E_part_gb + R_part_gb - L_part_gb + I_part_gb)*Xmat)
      
      ###########################################################
      ### Calculate pseudo hessian of beta at this iteration: ###
      ###########################################################
      num1 = exp(-c_haz)*(haz^2)
      num2 = exp(-c_haz)*haz
      num4 = exp(-c_haz_left)*haz_left
      
      comp1 = num1/exp1den 
      comp2 = num2/exp1den
      comp4 = num4/exp2den
      comp5 = num1/exp2den
      comp6 = num2/exp2den
      
      ### Fix numerical issues (by adding jitter to the denominator): ###
      exp4den = haz^2
      exp4den[exp4den == 0] = 1e-6
      
      # Pick out negative definite terms in original Hessian matrix and make them positive.
      # Compute each part in the pseudo hessian matrix for beta separately:
      E_part_hb = delE*((((haz_d^2)/exp4den)*(kp_right^2) + haz*kp_right)^delE)
      R_part_hb = delR*((haz*kp_right)^delR)
      L_part_hb = delL*(((comp1 + comp2^2)*(kp_right^2))^delL)
      I_part_hb = delI*(((comp4^2)*(kp_left_int^2) + comp4*kp_left_int + (comp5 + comp6^2)*(kp_right^2))^delI)
      
      pseudo_hess_beta = -t(Xmat)%*%(as.vector(E_part_hb + R_part_hb + L_part_hb + I_part_hb)*Xmat)
      
      ##########################################
      ### Backtracking line search for beta: ###
      ##########################################
      alpha_beta = 1
      factor = 0.5
      
      search = solve(pseudo_hess_beta)%*%grad_beta
      beta_new = as.vector(beta_old - alpha_beta*search)
      
      old_val = penloglikfun(beta_old, gVal, tVal, Xmat, rep_knots, rep_sd,
                             h_Rmat, id_Z, Z_t, time_int, mod_time_int, 
                             delE, delL, delR, delI, zfin)
      
      # Update initial values for beta if numerical issues occur:
      if(is.finite(old_val) == FALSE){
        cvg = 0
        msg = "Old penalized log-likelihood is not finite (line search)."
        oldVal_flag = TRUE
        print(msg)
        break
      }
      
      suppressWarnings({
        new_val = penloglikfun(beta_new, gVal, tVal, Xmat, rep_knots, rep_sd, 
                               h_Rmat, id_Z, Z_t, time_int, mod_time_int, 
                               delE, delL, delR, delI, zfin)
        
        while(is.finite(new_val) == FALSE){
          alpha_beta = alpha_beta*factor
          beta_new = as.vector(beta_old - alpha_beta*search)
          new_val = penloglikfun(beta_new, gVal, tVal, Xmat, rep_knots, rep_sd, 
                                 h_Rmat, id_Z, Z_t, time_int, mod_time_int,
                                 delE, delL, delR, delI, zfin)
        }
      })
      
      update_cond = (new_val >= old_val)
      
      if(update_cond == FALSE){
        alpha_beta = alpha_beta*factor
        beta_new = as.vector(beta_old - alpha_beta*search)
        cost = penloglikfun(beta_new, gVal, tVal, Xmat, rep_knots, rep_sd, h_Rmat,
                            id_Z, Z_t, time_int, mod_time_int, delE, delL, delR, delI, zfin)
        update_cond = (cost >= old_val)
      }
      
      # Store updated value of beta at this iteration: 
      bVal = beta_new
      ###########################################################################
      
      ########################## Newton step for gamma ##########################
      gamma_old = gVal
      gamma_new = {}
      
      # Calculate kappa values at this iteration: 
      kpset = kpvec_fun(bVal, gamma_old, Xmat, Z_t, id_Z, time_int, mod_time_int)
      
      kp_right = kpset$kp_right
      kp_d_right = kpset$kp_d_right
      kp_2d_right = kpset$kp_2d_right
      
      kp_left_int = kpset$kp_left_int
      kp_d_left_int = kpset$kp_d_left_int
      kp_2d_left_int = kpset$kp_2d_left_int
      
      # Baseline hazard and cumulative baseline hazard (most cases):
      set1 = gauss_basis_fun(tVal, kp_right, rep_knots, rep_sd)
      c_haz = set1$c_haz_gen
      haz = set1$haz_gen
      haz_d = set1$haz_d_gen
      
      # Cumulative baseline hazard (When considering the left point in int. censoring):
      set2 = gauss_basis_fun(tVal, kp_left_int, rep_knots, rep_sd)
      c_haz_left = set2$c_haz_gen
      haz_left = set2$haz_gen
      
      ### Fix numerical issues (by adding jitter to the denominator): ###
      
      exp1den = -expm1(-c_haz)  #exp1den = 1-exp(-c_haz)
      exp1den[exp1den == 0] = 1e-6
      
      exp2den = (expm1(-c_haz_left) - expm1(-c_haz))^delI  #exp2den = (exp(-c_haz_left) - exp(-c_haz))^delI
      exp2den[exp2den == 0] = 1e-6
      
      exp3den = haz
      exp3den[exp3den == 0] = 1e-6
      
      ######################################################
      ### Calculate gradient of gamma at this iteration: ###
      ######################################################
      expdLall = exp(-c_haz_left)*haz_left*kp_d_left_int
      expdRall = exp(-c_haz)*haz*kp_d_right
      
      # Compute each part in the gradient of gamma separately:
      E_part_gg = delE*(-zfin + (haz_d/exp3den - haz)*kp_d_right)
      R_part_gg = delR*(haz*kp_d_right)
      L_part_gg = delL*((exp(-c_haz)*haz)/exp1den)*kp_d_right
      I_part_gg = delI*(-expdLall + expdRall)/exp2den
      
      grad_gamma = colSums(E_part_gg - R_part_gg + L_part_gg + I_part_gg)
      
      ############################################################
      ### Calculate pseudo hessian of gamma at this iteration: ###
      ############################################################
      num1 = exp(-c_haz)*haz
      num2 = exp(-c_haz)*haz^2
      num3 = exp(-c_haz_left)*haz_left
      
      comp2 = num3/exp2den
      comp3 = num2/exp2den
      comp4 = num1/exp2den
      comp5 = (num1*num3)/(exp2den^2)
      
      fin1 = kp_d_left_int*(comp2^2)*kp_d_left_int
      fin2 = comp2*kp_2d_left_int
      fin3 = kp_d_right*(comp3 + comp4^2)*kp_d_right
      
      ### Fix numerical issues (by adding jitter to the denominator): ###
      exp4den = haz^2
      exp4den[exp4den == 0] = 1e-6
      
      # Pick out negative definite terms in original Hessian matrix and make them positive.
      # Compute each part in the pseudo hessian matrix for gamma separately:
      
      E_part_hg = delE*((kp_d_right*((haz_d^2)/exp4den)*kp_d_right + haz*kp_2d_right)^delE)
      
      R_part_hg = delR*((haz*kp_2d_right)^delR)
      
      L_part_hg = delL*((kp_d_right*((exp(-c_haz)*(haz^2))/exp1den 
                                     + ((exp(-c_haz)*haz)/exp1den)^2)*kp_d_right)^delL)
      
      I_part_hg = delI*((fin1 + fin2 + fin3)^delI)
      
      pseudo_hess_gamma = -colSums(E_part_hg + R_part_hg + L_part_hg + I_part_hg)
      
      ###########################################
      ### Backtracking line search for gamma: ###
      ###########################################
      alpha_gamma = 1
      factor = 0.5
      
      search = solve(pseudo_hess_gamma)%*%grad_gamma
      gamma_new = as.vector(gamma_old - alpha_gamma*search)
      
      old_val = penloglikfun(bVal, gamma_old, tVal, Xmat, rep_knots, rep_sd, h_Rmat,
                             id_Z, Z_t, time_int, mod_time_int,
                             delE, delL, delR, delI, zfin)
      
      suppressWarnings({
        new_val = penloglikfun(bVal, gamma_new, tVal, Xmat, rep_knots, rep_sd, h_Rmat,
                               id_Z, Z_t, time_int, mod_time_int, delE, delL, delR, delI, zfin)
        
        while(is.finite(new_val) == FALSE){
          alpha_gamma = alpha_gamma*factor
          gamma_new = as.vector(gamma_old - alpha_gamma*search)
          new_val = penloglikfun(bVal, gamma_new, tVal, Xmat, rep_knots, rep_sd, h_Rmat,
                                 id_Z, Z_t, time_int, mod_time_int, delE, delL, delR, delI, zfin)
        }
      })
      
      update_cond = (new_val >= old_val)
      
      while(update_cond == FALSE){
        alpha_gamma = alpha_gamma*factor
        gamma_new = as.vector(gamma_old - alpha_gamma*search)
        cost = penloglikfun(bVal, gamma_new, tVal, Xmat, rep_knots, rep_sd, h_Rmat,
                            id_Z, Z_t, time_int, mod_time_int, delE, delL, delR, delI, zfin)
        update_cond = (cost >= old_val)
      }
      
      ### Store updated value of gamma at this iteration: ###
      gVal = gamma_new
      ###########################################################################
      
      ####################### MI algorithm step for theta #######################
      theta_old = tVal
      theta_new = {}
      
      ### Calculate kappa values at this iteration: ###
      kpset = kpvec_fun(bVal, gVal, Xmat, Z_t, id_Z,  time_int, mod_time_int)
      
      kp_right = kpset$kp_right
      kp_d_right = kpset$kp_d_right
      kp_2d_right = kpset$kp_2d_right
      
      kp_left_int = kpset$kp_left_int
      kp_d_left_int = kpset$kp_d_left_int
      kp_2d_left_int = kpset$kp_2d_left_int
      
      ### Calculate hazard values at this iteration: ###
      # Baseline hazard and cumulative baseline hazard (most cases):
      set1 = gauss_basis_fun(theta_old, kp_right, rep_knots, rep_sd)
      c_haz = set1$c_haz_gen
      haz = set1$haz_gen
      haz_d = set1$haz_d_gen
      haz_2d = set1$haz_2d_gen
      gsb_all = t(set1$gsb)
      cumgsb_all = t(set1$cumgsb)
      
      # Cumulative baseline hazard (When considering the left point in int. censoring):
      set2 = gauss_basis_fun(theta_old, kp_left_int, rep_knots, rep_sd)
      c_haz_left = set2$c_haz_gen
      haz_left = set2$haz_gen
      gsb_left = t(set2$gsb)
      cumgsb_left = t(set2$cumgsb)
      
      ### Fix numerical issues (by adding jitter to the denominator): ###
      
      exp1den = -expm1(-c_haz)  #exp1den = 1-exp(-c_haz)
      exp1den[exp1den == 0] = 1e-6
      
      
      exp2den = (expm1(-c_haz_left) - expm1(-c_haz))^delI  #exp2den = (exp(-c_haz_left) - exp(-c_haz))^delI
      exp2den[exp2den == 0] = 1e-6
      
      exp3den = haz
      exp3den[exp3den == 0] = 1e-6
      
      ######################################################
      ### Calculate gradient of theta at this iteration: ###
      ######################################################
      
      # Compute each part in the gradient of theta separately:
      E_part_gt = delE*((gsb_all/exp3den - cumgsb_all)^delE)
      R_part_gt = delR*((cumgsb_all)^delR)
      L_part_gt = delL*(((cumgsb_all*exp(-c_haz))/exp1den)^delL)
      I_part_gt = delI*(((-cumgsb_left*exp(-c_haz_left) 
                          + cumgsb_all*exp(-c_haz))/exp2den)^delI)
      
      penalty_d_gt = 2*h_Rmat%*%theta_old
      
      grad_theta = colSums(E_part_gt - R_part_gt + L_part_gt + I_part_gt) - penalty_d_gt
      
      #############################################
      ### Calculate MI ratio at this iteration: ###
      #############################################
      
      penalty_pos = h_Rmat%*%theta_old
      penalty_pos[penalty_pos < 0] = 0
      
      S_theta_den = as.vector(colSums(delE*cumgsb_all
                                      + delR*cumgsb_all
                                      + delI*((cumgsb_left*exp(-c_haz_left))/exp2den))
                              + 2*penalty_pos)
      
      S_theta = diag(theta_old/S_theta_den)
      
      ### Backtracking line search for theta: ###
      
      alpha_theta = 1
      factor = 0.5
      
      search = S_theta%*%grad_theta
      theta_new = as.vector(theta_old + alpha_theta*search)
      
      old_val = penloglikfun(bVal, gVal, theta_old, Xmat, rep_knots, rep_sd, h_Rmat,
                             id_Z, Z_t, time_int, mod_time_int, delE, delL, delR, delI, zfin)
      
      suppressWarnings({
        new_val = penloglikfun(bVal, gVal, theta_new, Xmat, rep_knots, rep_sd, h_Rmat,
                               id_Z, Z_t, time_int, mod_time_int, delE, delL, delR, delI, zfin)
        
        while(is.finite(new_val) == FALSE){
          alpha_theta = alpha_theta*factor
          theta_new = as.vector(theta_old + alpha_theta*search)
          new_val = penloglikfun(bVal, gVal, theta_new, Xmat, rep_knots, rep_sd, h_Rmat,
                                 id_Z, Z_t, time_int, mod_time_int, delE, delL, delR, delI, zfin)
        }
      })
      
      update_cond = (new_val >= old_val)
      
      while(update_cond == FALSE){
        alpha_theta = alpha_theta*factor
        theta_new = as.vector(theta_old + alpha_theta*search)
        cost = penloglikfun(bVal, gVal, theta_new, Xmat, rep_knots, rep_sd, h_Rmat,
                            id_Z, Z_t, time_int, mod_time_int, delE, delL, delR, delI, zfin)
        update_cond = (cost >= old_val)
      }
      
      ### Store updated value of theta at this iteration: ###
      
      tVal = theta_new
      
      ###########################################################################
      
      ######################### Convergence criteria checks: ####################
      
      ### Calculate and store value of penalized log-likelihood: ###
      penloglik = penloglikfun(bVal, gVal, tVal, Xmat, rep_knots, rep_sd, h_Rmat, id_Z,
                               Z_t, time_int, mod_time_int, delE, delL, delR, delI, zfin)
      
      logLikVec = c(logLikVec, penloglik)
      
      ### Computing the gradient values after all updates: ###
      grad = grad_fun(bVal, gVal, tVal, rep_knots, rep_sd, h_Rmat, Xmat,
                      id_Z, Z_t, time_int, mod_time_int, delE, delL, delR, delI, zfin)
      
      # Code for convergence metric:
      # Check that new values of parameters and old values are within a pre-set tolerance level
      if(all(abs(beta_new-beta_old) < tol1) 
         & all(abs(gamma_new-gamma_old) < tol1) 
         & all(abs(theta_new-theta_old) < tol1)
         & all(abs(grad[1:3]) < tol2))
        #& all(abs(grad[4: (m + 3)][theta_new > 0.01]) < tol2))
      {break}
      
      #if((abs(logLikVec[iter+1] - logLikVec[iter]) < tol1)){
      #  logLikCount = logLikCount + 1
      #}
      #else{logLikCount = 0}
      
      #if((logLikCount >= 5)
      #   & all(abs(theta_new-theta_old) < tol2)){break}
      
      # Stopping condition to prevent an infinite loop:
      if(iterLoop >= maxIterPerLoop)
      {msg = "This inner loop did not converge. Move on to outer loop."
      print(msg)
      break}
      
      if(iter %% 100 == 0){print(iter)}
      ###########################################################################
    }
    
    ############################## End of inner loop: #########################
    #print(iter)
    
    # Break out of outer loop if numerical issues are present:
    if(oldVal_flag == TRUE){break}
    
    ######################### Begin automatic smoothing ######################### 
    
    # Code for updating the smoothing parameters (h):
    
    if(smooth_stop == FALSE){
      
      Qmat = matrix(rep(0, numPar^2), numPar, numPar)
      Qmat[(bgDim + 1):numPar, (bgDim + 1):numPar] = 2*h_Rmat
      
      Fmat = -fullHessian(bVal, gVal, tVal, rep_knots, rep_sd, Xmat, Z_t,
                          id_Z, time_int, mod_time_int,delE, delL, delR, delI)
      
      nu = sum(diag(solve(Fmat + Qmat)%*%Qmat))
      nu_record = rbind(nu_record, nu)
      
      # Re-set stabilty count if nu is negative:
      if(nu < 0){nuStableCount = 0}
      
      # Smoothing parameter updates and converges if stability count reaches 3:
      if((abs(nu - nu_old) >= diffDF)){
        
        nuStableCount = 0
        nu_old = nu
        
        df_diff = m - nu
        if(df_diff <= 0){df_diff = 1}
        
        sigma2_h = (tVal%*%Rmat%*%tVal)/df_diff
        h = as.numeric(1/(2*sigma2_h))
        h_record = rbind(h_record, h)
        
        count = count + 1
      }
      else if((abs(nu - nu_old) < diffDF) & nuStableCount < stableNum){
        
        nuStableCount = nuStableCount + 1
        nu_old = nu
        
        df_diff = m - nu
        if(df_diff <= 0){df_diff = 1}
        
        sigma2_h = (tVal%*%Rmat%*%tVal)/df_diff
        h = as.numeric(1/(2*sigma2_h))
        h_record = rbind(h_record, h)
        
        count = count + 1
      }
      else{
        cvg = 1
        h_record = rbind(h_record, NA)
        msg = "This iteration converged! (smoothing parameters updated successfully)"
        print(msg)
        break
      }
      
      # Break code if more than 10 outer loops are required:
      if(count == outerMax + 1){
        cvg = 0
        msg = "This replication exceeded maximum number outer loops."
        print(msg)
        break}
    }
    else{cvg = 1
    msg = "This iteration converged! (No automatic smoothing implemented)"
    print(msg)
    break}
    
    ######################### End of automatic smoothing ######################### 
  }
  
  ############################## End of outer loop: ############################
  
  ###################### Calculating the asymptotic variance ###################
  
  Fmat =  -fullHessian(bVal, gVal, tVal, rep_knots, rep_sd, Xmat, Z_t,
                       id_Z, time_int, mod_time_int, delE, delL, delR, delI)
  
  if(any(is.na(Fmat)) == TRUE){
    cvg = 0
    msg = "This iteration converged but unable to evaluate Hessian"
    print(msg)
  }
  
  thetaGrad = grad[-c(1:bgDim)]
  activeTheta = ((tVal < 0.01) & (thetaGrad < 0))*1
  
  Qmat = matrix(rep(0, numPar^2), numPar, numPar)
  Qmat[(bgDim + 1):numPar, (bgDim + 1):numPar] = 2*h*Rmat
  
  # Identify which rows (active constraints) need to be cut:
  if(sum(activeTheta) > 0){
    startMat = diag(numPar)
    removeCols = which(activeTheta == 1) + bgDim
    Umat = startMat[,-removeCols]
  }
  else{Umat = diag(numPar)}
  
  # Code to reset negative eigenvalues to small positive eigenvalues:
  eigenVals = eigen(Fmat)$values
  
  if(any(eigenVals < 0) == TRUE){
    
    eigenVals[eigenVals < 0] = runif(sum(eigenVals < 0), min = 1.1E-6, max = 1.5E-6)
    eigenVecs = eigen(Fmat)$vectors
    
    FmatEigen = eigenVecs %*% diag(eigenVals) %*% t(eigenVecs)
    
    BMatEigenInv = Umat%*%solve(t(Umat)%*%(FmatEigen + Qmat)%*%Umat)%*%t(Umat)
    
    asyEigenVar = BMatEigenInv%*%FmatEigen%*%BMatEigenInv
    
  }
  else{
    BMatInv = Umat%*%solve(t(Umat)%*%(Fmat + Qmat)%*%Umat)%*%t(Umat)
    
    asyEigenVar = BMatInv%*%Fmat%*%BMatInv
  }
  
  return(list(bVal = bVal, gVal = gVal, tVal = tVal, grad = grad, h = h,
              iter = iter, count = count, rep_knots = rep_knots, rep_sd = rep_sd,
              haz = haz, asyEigenVar = asyEigenVar, cvg = cvg, msg = msg, 
              logLikVec = logLikVec, nu_record = nu_record, h_record = h_record, 
              kappa_vec = kappa_vec))
}
