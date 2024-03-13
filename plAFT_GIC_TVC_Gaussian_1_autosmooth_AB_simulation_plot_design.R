# forked from plAFT_GIC_TVC_Gaussian_1_autosmooth_WBRT_plot_design.R on 20240130
plAFT_GIC_TVC_Gaussian_1_autosmooth_AB_simulation_plot_design <- function(data, estimation_result, resolution = 500, knots_option = "percentile", plot_sig_lv = 0.05) {
  # browser()
  # tmat = as.matrix(data$tmat)
  # rownames(tmat) = NULL
  # time = c(tmat[,1],tmat[,2])
  # time = time[!is.infinite(time)]
  # range_time = range(time)
  range_time = c(0, data$max_time)
  # range_time = c(0, 1.25)
  # censor_n = tmat[, 3, drop = FALSE] # n by 1 matrix
  censor_n = data$max_time
  
  Xmat = as.matrix(data$Xmat)
  rownames(Xmat) = NULL
  X = Xmat[, -1, drop = FALSE]
  n = dim(X)[1]
  p = dim(X)[2]
  
  
  beta_hat = estimation_result$beta_hat
  gamma_hat = estimation_result$gamma_hat
  theta_hat = estimation_result$theta_hat
  asym_std_beta = estimation_result$asym_std_beta
  asym_std_gamma = estimation_result$asym_std_gamma
  asym_std_theta = estimation_result$asym_std_theta
  
  
  # gknots_plot = estimation_result$knots_iter[NROW(estimation_result$knots_iter),]
  gknots_plot = estimation_result$knots_location
  sig_plot = estimation_result$knots_sigma
  cov_mat_theta = estimation_result$cov_mat_theta
  index_active_theta = estimation_result$index_active_theta
  
  # Gaussian_basis <- function(x, mean, sd) {
  #   if (nrow(x)==1) {return(matrix(sapply(X = 1:length(mean), FUN = function(i) {dnorm(x, mean[i], sd[i]) * sd[i]*sqrt(2*pi)}), nrow = 1))}
  #   else {return(sapply(X = 1:length(mean), FUN = function(i) {dnorm(x, mean[i], sd[i]) * sd[i]*sqrt(2*pi)}))}
  # }
  # Gaussian_basis_integ1 <- function(x, mean, sd) {
  #   if (nrow(x)==1) {return(matrix(sapply(X = 1:length(mean), FUN = function(i) {(pnorm(x, mean[i], sd[i]) - matrix(1, NROW(x), 1)%*%pnorm(0, mean[i], sd[i])) * sd[i]*sqrt(2*pi)}), nrow = 1))}
  #   else {return(sapply(X = 1:length(mean), FUN = function(i) {(pnorm(x, mean[i], sd[i]) - matrix(1, NROW(x), 1)%*%pnorm(0, mean[i], sd[i])) * sd[i]*sqrt(2*pi)}))}
  # }
  
  # AB's way, I modified on 20240201
  Gaussian_basis <- function(x, mean, sd) {
    if (nrow(x)==1) {return(matrix(sapply(X = 1:length(mean), FUN = function(i) {dnorm(x, mean[i], sd[i])}), nrow = 1))}
    else {return(sapply(X = 1:length(mean), FUN = function(i) {dnorm(x, mean[i], sd[i])}))}
  }
  Gaussian_basis_integ1 <- function(x, mean, sd) {
    if (nrow(x)==1) {return(matrix(sapply(X = 1:length(mean), FUN = function(i) {(pnorm(x, mean[i], sd[i]) - matrix(1, NROW(x), 1)%*%pnorm(0, mean[i], sd[i]))}), nrow = 1))}
    else {return(sapply(X = 1:length(mean), FUN = function(i) {(pnorm(x, mean[i], sd[i]) - matrix(1, NROW(x), 1)%*%pnorm(0, mean[i], sd[i]))}))}
  }
  
  # calculation based on X and Z_5
  
  time_plot = seq(1e-10, range_time[2], length.out = resolution)
  
  # create the Zmat based on the records of {0, 0}
  Zmat_00 = NULL
  for (i in 1:length(time_plot)) {
    if (time_plot[i] < 0.2) {
      Zmat_00 = rbind(Zmat_00,
                      c(i, 1, 0, time_plot[i], 0))
    }
    else {
      Zmat_00 = rbind(Zmat_00,
                      c(i, 1, 0, 0.2, 0),
                      c(i, 2, 0.2, time_plot[i], 0))
    }
  }
  rownames(Zmat_00) = NULL
  id_Z_00 = Zmat_00[,1]
  interval_length_00 = Zmat_00[, 4, drop = FALSE] - Zmat_00[, 3, drop = FALSE] # N by 1 matrix
  Z_00 = Zmat_00[, -c(1:4), drop = FALSE] # N by q matrix
  Z_00_ni = Z_00[c(diff(id_Z_00)!=0,TRUE), , drop = FALSE]
  
  # create the Zmat based on the records of {0, 1}
  Zmat_01 = NULL
  for (i in 1:length(time_plot)) {
    if (time_plot[i] < 0.2) {
      Zmat_01 = rbind(Zmat_01,
                      c(i, 1, 0, time_plot[i], 0))
    }
    else {
      Zmat_01 = rbind(Zmat_01,
                      c(i, 1, 0, 0.2, 0),
                      c(i, 2, 0.2, time_plot[i], 1))
    }
  }
  rownames(Zmat_01) = NULL
  id_Z_01 = Zmat_01[,1]
  interval_length_01 = Zmat_01[, 4, drop = FALSE] - Zmat_01[, 3, drop = FALSE] # N by 1 matrix
  Z_01 = Zmat_01[, -c(1:4), drop = FALSE] # N by q matrix
  Z_01_ni = Z_01[c(diff(id_Z_01)!=0,TRUE), , drop = FALSE]
  
  
  plot_ingredients <- function(X_input, beta = beta_hat, Z_input, Z_ni_input, gamma = gamma_hat, theta = theta_hat, interval_length_input, id_Z_input, cov_mat_theta_input = cov_mat_theta, index_active_theta_input = index_active_theta) {
    # browser()
    ts_plot = exp(-X_input%*%beta) * as.matrix(tapply(exp(-Z_input%*%gamma)*interval_length_input, id_Z_input, FUN = sum))
    
    
    # ts = exp(-X_input%*%beta) * as.matrix(tapply(exp(-Z_input%*%gamma)*interval_length_input, id_Z_input, FUN = sum))
    if (knots_option=="equal_space") {
      # bryKnt = range(ts_plot)
      bryKnt = range(ts_plot)
      bin_width = (bryKnt[2] - bryKnt[1]) / (length(theta_hat)-1)
      knots = seq(bryKnt[1], bryKnt[2], bin_width)
      sigs = rep((2/3) * bin_width, times = length(theta_hat))
    }
    else {
      knots = gknots_plot
      sigs = sig_plot
    }
    
    # ts_plot = matrix(seq(1e-10, max(ts), length.out = resolution))
    phi_plot = Gaussian_basis(x = ts_plot, mean = knots, sd = sigs)
    Phi_plot = Gaussian_basis_integ1(x = ts_plot, mean = knots, sd = sigs)
    
    if (length(index_active_theta_input)>0) {
      cov_mat_theta_no_active = cov_mat_theta_input[-index_active_theta_input,-index_active_theta_input]
      phi_plot_no_active = phi_plot[,-index_active_theta_input]
      h0ts_plot_no_active = phi_plot_no_active%*%theta[-index_active_theta_input]
      Phi_plot_no_active = Phi_plot[,-index_active_theta_input]
      S0ts_plot_no_active = exp(-Phi_plot_no_active%*%theta[-index_active_theta_input])
    }
    else {
      cov_mat_theta_no_active = cov_mat_theta_input
      phi_plot_no_active = phi_plot
      h0ts_plot_no_active = phi_plot%*%theta
      Phi_plot_no_active = Phi_plot
      S0ts_plot_no_active = exp(-Phi_plot%*%theta)
    }
    
    # baseline hazard plot of accelerated/decelerated time scale
    G_h0ts_plot = log(h0ts_plot_no_active)
    dG_dtheta_h0ts_plot = (1/h0ts_plot_no_active)%*%matrix(1,1,ncol(phi_plot_no_active)) * phi_plot_no_active
    asym_std_G_h0ts_plot = sqrt(diag(dG_dtheta_h0ts_plot%*%cov_mat_theta_no_active%*%t(dG_dtheta_h0ts_plot)))
    CI_G_h0ts_plot = cbind(G_h0ts_plot-qnorm(1-plot_sig_lv/2)*asym_std_G_h0ts_plot, G_h0ts_plot+qnorm(1-plot_sig_lv/2)*asym_std_G_h0ts_plot)
    CI_h0ts_plot = exp(CI_G_h0ts_plot)
    
    ht_plot = h0ts_plot_no_active * exp(-X_input%*%beta - Z_ni_input%*%gamma)
    G_ht_plot = log(ht_plot)
    asym_std_G_ht_plot = asym_std_G_h0ts_plot
    CI_G_ht_plot = cbind(G_ht_plot-qnorm(1-plot_sig_lv/2)*asym_std_G_ht_plot, G_ht_plot+qnorm(1-plot_sig_lv/2)*asym_std_G_ht_plot)
    CI_ht_plot = exp(CI_G_ht_plot)
    
    # baseline survival plot of accelerated/decelerated time scale
    G_S0ts_plot = log(S0ts_plot_no_active/(1-S0ts_plot_no_active))
    dG_dtheta_S0ts_plot = -(1/(1-S0ts_plot_no_active))%*%matrix(1,1,ncol(Phi_plot_no_active)) * Phi_plot_no_active
    asym_std_G_S0ts_plot = sqrt(diag(dG_dtheta_S0ts_plot%*%cov_mat_theta_no_active%*%t(dG_dtheta_S0ts_plot)))
    CI_G_S0ts_plot = cbind(G_S0ts_plot-qnorm(1-plot_sig_lv/2)*asym_std_G_S0ts_plot, G_S0ts_plot+qnorm(1-plot_sig_lv/2)*asym_std_G_S0ts_plot)
    CI_S0ts_plot = 1/(exp(-CI_G_S0ts_plot)+1)
    
    results_plot_ingredients = list()
    results_plot_ingredients$ts_plot = ts_plot
    
    results_plot_ingredients$h0ts_plot = h0ts_plot_no_active
    results_plot_ingredients$CI_h0ts_plot = CI_h0ts_plot
    
    results_plot_ingredients$ht_plot = ht_plot
    results_plot_ingredients$CI_ht_plot = CI_ht_plot
    
    results_plot_ingredients$S0ts_plot = S0ts_plot_no_active
    results_plot_ingredients$CI_S0ts_plot = CI_S0ts_plot
    
    return(results_plot_ingredients)
  }
  
  # design X input
  # browser()
  {
    X_x1_1 = matrix(1, resolution, 2)
    X_x1_1[, 2] = median(X[, 2])
    option_x1_1_z_00 = plot_ingredients(X_input = X_x1_1, Z_input = Z_00, Z_ni_input = Z_00_ni, interval_length_input = interval_length_00, id_Z_input = id_Z_00)
    option_x1_1_z_01 = plot_ingredients(X_input = X_x1_1, Z_input = Z_01, Z_ni_input = Z_01_ni, interval_length_input = interval_length_01, id_Z_input = id_Z_01)
  }
  
  {
    X_x1_0 = matrix(0, resolution, 2)
    X_x1_0[, 2] = median(X[, 2])
    option_x1_0_z_00 = plot_ingredients(X_input = X_x1_0, Z_input = Z_00, Z_ni_input = Z_00_ni, interval_length_input = interval_length_00, id_Z_input = id_Z_00)
    option_x1_0_z_01 = plot_ingredients(X_input = X_x1_0, Z_input = Z_01, Z_ni_input = Z_01_ni, interval_length_input = interval_length_01, id_Z_input = id_Z_01)
  }
  
  plot(time_plot, option_x1_1_z_00$S0ts_plot, xlim = c(0,range_time[2]*1.05), ylim = c(0,1), xlab = expression(t), ylab = expression(S(t)), mgp=c(2,1,0), type = "l", lty = "solid", col = "white")
  grid(lty = "solid")
  lines(time_plot, option_x1_1_z_00$S0ts_plot, type = "l", lty = "solid", lwd = 2, col = "black")
  lines(time_plot, option_x1_0_z_00$S0ts_plot, type = "l", lty = "solid", lwd = 2, col = "darkgrey")

  lines(time_plot, option_x1_1_z_01$S0ts_plot, type = "l", lty = "dashed", lwd = 2, col = "black")
  lines(time_plot, option_x1_0_z_01$S0ts_plot, type = "l", lty = "dashed", lwd = 2, col = "darkgrey")

  abline(v = 0.2, lty = "dashed", lwd = 2, col = "red")

  title(main = "Predicted Survival (x1)")
  legend("topright", legend = c("x1_1_z_00", "x1_0_z_00", "x1_1_z_01", "x1_0_z_01"), col = c("black","darkgrey","black","darkgrey"), lty = c("solid","solid","dashed","dashed"), lwd = 2, cex = 0.85)
  
  # x1_1_z_01 vs x1_0_z_01
  # survival plot
  plot(time_plot, option_x1_1_z_01$S0ts_plot, xlim = c(0,range_time[2]*1.05), ylim = c(0,1), xlab = expression(t), ylab = expression(S(t)), mgp=c(2,1,0), type = "l", lty = "solid", col = "white")
  grid(lty = "solid")
  lines(time_plot, option_x1_1_z_01$S0ts_plot, type = "l", lty = "solid", lwd = 2, col = "black")
  lines(time_plot, option_x1_0_z_01$S0ts_plot, type = "l", lty = "dotdash", lwd = 2, col = "black")
  abline(v = 0.2, lty = "dashed", lwd = 2, col = "red")
  text(0.12, 0.1, expression(tau==0.2), srt = 90, col = "red", cex = 0.8)
  title(main = "Predicted Survival")
  legend("topright", legend = c(expression(S(t~"|"~x[1]==1,x[2],tilde(bold(z))(t),tau)), expression(S(t~"|"~x[1]==0,x[2],tilde(bold(z))(t),tau))), col = c("black","black"), lty = c("solid","dotdash"), text.width = 1.2, lwd = 2, cex = 1)
  
  # survival ratio plot
  par(mar = c(5.1, 5.1, 4.1, 2.1))
  # plot(time_plot, option_x1_1_z_01$S0ts_plot/option_x1_0_z_01$S0ts_plot, xlim = c(0,range_time[2]*1.05), xlab = expression(t), ylab = expression(frac(S(t~"|"~x[1]==1,x[2],tilde(bold(z))(t),tau),S(t~"|"~x[1]==0,x[2],tilde(bold(z))(t),tau))), mgp=c(2,1,0), type = "l", lty = "solid", col = "black")
  plot(time_plot, option_x1_0_z_01$S0ts_plot/option_x1_1_z_01$S0ts_plot, xlim = c(0,range_time[2]*1.05), xlab = expression(t), ylab = expression(frac(S(t~"|"~x[1]==0,x[2],tilde(bold(z))(t),tau),S(t~"|"~x[1]==1,x[2],tilde(bold(z))(t),tau))), mgp=c(2,1,0), type = "l", lty = "solid", col = "white")
  # plot(time_plot, option_x1_1_z_01$S0ts_plot/option_x1_0_z_01$S0ts_plot, xlim = c(0,range_time[2]*1.05), xlab = "Survival Ratio", mgp=c(2,1,0), type = "l", lty = "solid", col = "black")
  # plot(time_plot, option_x1_1_z_01$S0ts_plot/option_x1_0_z_01$S0ts_plot, xlim = c(0,0.3), ylim = c(0, 3), xlab = expression(t), ylab = expression(S(t)), mgp=c(2,1,0), type = "l", lty = "solid", col = "black")
  grid(lty = "solid")
  lines(time_plot, option_x1_0_z_01$S0ts_plot/option_x1_1_z_01$S0ts_plot, type = "l", lty = "solid", lwd = 2, col = "black")
  abline(h = 1, lty = "dashed", lwd = 2, col = "black")
  abline(v = 0.2, lty = "dashed", lwd = 2, col = "red")
  text(0.12, 0.1, expression(tau==0.2), srt = 90, col = "red", cex = 0.8)
  title(main = "Predicted Survival Ratio")
  # title(main = expression("Predicted"~"Survival"~"Ratio"~frac(S(t~"|"~x[1]==1,x[2],tilde(bold(z))(t),tau), S(t~"|"~x[1]==0,x[2],tilde(bold(z))(t),tau))))
  
  # x1_1_z_01 vs x1_1_z_00
  # survival plot
  par(mar = c(5.1, 4.1, 4.1, 2.1))
  plot(time_plot, option_x1_1_z_01$S0ts_plot, xlim = c(0,range_time[2]*1.05), ylim = c(0,1), xlab = expression(t), ylab = expression(S(t)), mgp=c(2,1,0), type = "l", lty = "solid", col = "white")
  grid(lty = "solid")
  lines(time_plot, option_x1_1_z_01$S0ts_plot, type = "l", lty = "solid", lwd = 2, col = "black")
  lines(time_plot, option_x1_1_z_00$S0ts_plot, type = "l", lty = "dotdash", lwd = 2, col = "black")
  abline(v = 0.2, lty = "dashed", lwd = 2, col = "red")
  text(0.12, 0.1, expression(tau==0.2), srt = 90, col = "red", cex = 0.8)
  title(main = "Predicted Survival")
  legend("topright", legend = c(expression(S(t~"|"~x[1]==1,x[2],tilde(bold(z))(t),tau)), expression(S(t~"|"~x[1]==1,x[2],tilde(bold(z))(t)==bold(0)))), col = c("black","black"), lty = c("solid","dotdash"), lwd = 2, text.width = 1.2, cex = 1)
  # survival ratio plot
  # plot(time_plot, option_x1_1_z_01$S0ts_plot/option_x1_1_z_00$S0ts_plot, xlim = c(0,range_time[2]*1.05), xlab = expression(t), ylab = expression(S(t~"|"~x[1]==1,x[2],tilde(bold(z))(t),tau)/S(t~"|"~x[1]==1,x[2])), mgp=c(2,1,0), type = "l", lty = "solid", col = "black")
  par(mar = c(5.1, 5.1, 4.1, 2.1))
  plot(time_plot, option_x1_1_z_01$S0ts_plot/option_x1_1_z_00$S0ts_plot, xlim = c(0,range_time[2]*1.05), xlab = expression(t), ylab = expression(frac(S(t~"|"~x[1]==1,x[2],tilde(bold(z))(t),tau),S(t~"|"~x[1]==1,x[2],tilde(bold(z))(t)==bold(0)))), mgp=c(2,1,0), type = "l", lty = "solid", col = "white")
  grid(lty = "solid")
  lines(time_plot, option_x1_1_z_01$S0ts_plot/option_x1_1_z_00$S0ts_plot, type = "l", lty = "solid", lwd = 2, col = "black")
  abline(h = 1, lty = "dashed", lwd = 2, col = "black")
  abline(v = 0.2, lty = "dashed", lwd = 2, col = "red")
  text(0.12, 0.1, expression(tau==0.2), srt = 90, col = "red", cex = 0.8)
  title(main = "Predicted Survival Ratio")
  par(mar = c(5.1, 4.1, 4.1, 2.1))
  # plot(time_plot, option_M_A_S_slow$h0ts_plot, xlim = c(0,max(time_plot,option_M_A_S_slow$ts_plot,option_M_O_S_slow$ts_plot)), ylim = c(0,max(option_M_A_S_slow$h0ts_plot,option_M_O_S_slow$h0ts_plot)), main = "hazard plot", xlab = expression(t), ylab = expression(lambda[0](kappa(t))), mgp=c(2,1,0), type = "l", lty = "solid", col = "white") # variables for the paper
  # # plot(option_hat$ts_plot, option_hat$h0ts_plot, xlim = c(0,max(option_hat$ts_plot,option_M_A_S_5$ts_plot,option_M_O_S_5$ts_plot)), ylim = c(0,max(option_hat$h0ts_plot,option_M_A_S_5$h0ts_plot,option_M_O_S_5$h0ts_plot)), main = "baseline hazard plot", xlab = expression(t), ylab = expression(lambda[0](t)), mgp=c(2,1,0), type = "l", lty = "solid") # variables for the paper
  # grid(lty = "solid")
  # # polygon(c(option_hat$ts_plot,rev(option_hat$ts_plot)), c(option_hat$CI_h0ts_plot[,1], rev(option_hat$CI_h0ts_plot[,2])), col = "#6BD7AF", lty = 0)
  # # lines(option_hat$ts_plot, option_hat$h0ts_plot, type = "l", lty = "solid", col = "black")
  # # lines(option_hat$ts_plot, option_hat$CI_h0ts_plot[,1], type = "l", lty = "dashed", col = "black")
  # # lines(option_hat$ts_plot, option_hat$CI_h0ts_plot[,2], type = "l", lty = "dashed", col = "black")
  # # legend("bottomright", legend = c("baseline hazard", paste0((1-plot_sig_lv)*100,"% ASYM CI")), col = c("black","black"), lty = c("solid","dashed"), cex = 0.6)
  # 
  # # lines(option_0$ts_plot, option_0$h0ts_plot, type = "l", lty = "dashed", col = "black")
  # # lines(option_M_A_S_5$ts_plot, option_M_A_S_5$h0ts_plot, type = "l", lty = "solid", col = "blue")
  # 
  # lines(time_plot, option_M_A_S_slow$h0ts_plot, type = "l", lty = "solid", col = "black")
  # lines(time_plot, option_M_O_S_slow$h0ts_plot, type = "l", lty = "dashed", col = "black")
  # lines(time_plot, option_M_A_S_normal$h0ts_plot, type = "l", lty = "solid", col = "blue")
  # lines(time_plot, option_M_O_S_normal$h0ts_plot, type = "l", lty = "dashed", col = "blue")
  # lines(time_plot, option_M_A_S_fast$h0ts_plot, type = "l", lty = "solid", col = "red")
  # lines(time_plot, option_M_O_S_fast$h0ts_plot, type = "l", lty = "dashed", col = "red")
  # 
  # 
  # plot(time_plot, option_M_A_S_slow$ht_plot, xlim = c(0,max(time_plot,option_M_A_S_slow$ts_plot,option_M_O_S_slow$ts_plot)), ylim = c(0,max(option_M_A_S_slow$ht_plot,option_M_O_S_slow$ht_plot)), main = "hazard plot", xlab = expression(t), ylab = expression(lambda(t)), mgp=c(2,1,0), type = "l", lty = "solid", col = "white") # variables for the paper
  # grid(lty = "solid")
  # # lines(time_plot, option_M_A_S_5$ht_plot, type = "l", lty = "solid", col = "black") # ALS
  # # lines(time_plot, option_M_A_S_5$CI_ht_plot[,1], type = "l", lty = "dashed", col = "black")
  # # lines(time_plot, option_M_A_S_5$CI_ht_plot[,2], type = "l", lty = "dashed", col = "black")
  # # lines(time_plot, option_M_O_S_5$ht_plot, type = "l", lty = "dashed", col = "black") # Others
  # # lines(time_plot, option_M_O_S_5$CI_ht_plot[,1], type = "l", lty = "dotdash", col = "blue")
  # # lines(time_plot, option_M_O_S_5$CI_ht_plot[,2], type = "l", lty = "dotdash", col = "blue")
  # # lines(time_plot, option_0$ht_plot, type = "l", lty = "dotdash", col = "red") # baseline?
  # # legend("topleft", legend = c("C_ALS", "Others"), col = c("black","blue"), lty = c("solid","dashed"), cex = 0.6)
  # # legend("topleft", legend = c("C_ALS", paste0((1-plot_sig_lv)*100,"% ASYM CI"), "Others", paste0((1-plot_sig_lv)*100,"% ASYM CI")), col = c("black","black","blue","blue"), lty = c("solid","dashed","dashed","dotdash"), cex = 0.6)
  # 
  # lines(time_plot, option_M_A_S_slow$ht_plot, type = "l", lty = "solid", col = "black")
  # lines(time_plot, option_M_O_S_slow$ht_plot, type = "l", lty = "dashed", col = "black")
  # lines(time_plot, option_M_A_S_normal$ht_plot, type = "l", lty = "solid", col = "blue")
  # lines(time_plot, option_M_O_S_normal$ht_plot, type = "l", lty = "dashed", col = "blue")
  # lines(time_plot, option_M_A_S_fast$ht_plot, type = "l", lty = "solid", col = "red")
  # lines(time_plot, option_M_O_S_fast$ht_plot, type = "l", lty = "dashed", col = "red")
  # legend("topleft", legend = c("SLOW C_ALS", "SLOW Others", "NORMAL C_ALS", "NORMAL Others", "FAST C_ALS", "FAST Others"), col = c("black","black","blue","blue","red","red"), lty = c("solid","dashed","solid","dashed","solid","dashed"), cex = 0.6)
  # 
  # 
  # plot(time_plot, option_M_A_S_slow$S0ts_plot, xlim = c(0,max(time_plot,option_M_A_S_slow$ts_plot,option_M_O_S_slow$ts_plot)), ylim = c(0,1), main = "survival plot", xlab = expression(t), ylab = expression(S(t)), mgp=c(2,1,0), type = "l", lty = "solid", col = "white") # variables for the paper
  # # plot(option_hat$ts_plot, option_hat$S0ts_plot, xlim = c(0,max(option_hat$ts_plot,option_M_A_S_5$ts_plot,option_M_O_S_5$ts_plot)), ylim = c(0,1), main = "baseline survival plot", xlab = expression(t), ylab = expression(S[0](t)), mgp=c(2,1,0), type = "l", lty = "solid") # variables for the paper
  # grid(lty = "solid")
  # # lines(time_plot, option_M_A_S_5$S0ts_plot, type = "l", lty = "solid", col = "black")
  # # lines(time_plot, option_M_O_S_5$S0ts_plot, type = "l", lty = "dashed", col = "black")
  # 
  # lines(time_plot, option_M_A_S_slow$S0ts_plot, type = "l", lty = "solid", col = "black")
  # lines(time_plot, option_M_O_S_slow$S0ts_plot, type = "l", lty = "dashed", col = "black")
  # lines(time_plot, option_M_A_S_normal$S0ts_plot, type = "l", lty = "solid", col = "blue")
  # lines(time_plot, option_M_O_S_normal$S0ts_plot, type = "l", lty = "dashed", col = "blue")
  # lines(time_plot, option_M_A_S_fast$S0ts_plot, type = "l", lty = "solid", col = "red")
  # lines(time_plot, option_M_O_S_fast$S0ts_plot, type = "l", lty = "dashed", col = "red")
  # legend("bottomleft", legend = c("S_ALSFRS_ALS", "S_ALSFRS_Others", "N_ALSFRS_ALS", "N_ALSFRS_Others", "F_ALSFRS_ALS", "F_ALSFRS_Others"), col = c("black","black","blue","blue","red","red"), lty = c("solid","dashed","solid","dashed","solid","dashed"), cex = 0.6)
  # 
  # # alternative colour theme
  # plot(time_plot, option_M_A_S_slow$S0ts_plot, xlim = c(0,max(time_plot,option_M_A_S_slow$ts_plot,option_M_O_S_slow$ts_plot)), ylim = c(0,1), main = "Predicted survival plot", xlab = expression(t), ylab = expression(S(t)), mgp=c(2,1,0), type = "l", lty = "solid", col = "white") # variables for the paper
  # grid(lty = "solid")
  # lines(time_plot, option_M_A_S_slow$S0ts_plot, type = "l", lty = "solid", col = "black")
  # lines(time_plot, option_M_O_S_slow$S0ts_plot, type = "l", lty = "solid", col = "darkgrey")
  # lines(time_plot, option_M_A_S_normal$S0ts_plot, type = "l", lty = "dashed", col = "black")
  # lines(time_plot, option_M_O_S_normal$S0ts_plot, type = "l", lty = "dashed", col = "darkgrey")
  # lines(time_plot, option_M_A_S_fast$S0ts_plot, type = "l", lty = "dotdash", col = "black")
  # lines(time_plot, option_M_O_S_fast$S0ts_plot, type = "l", lty = "dotdash", col = "darkgrey")
  # legend("bottomleft", legend = c("S_ALSFRS_ALS", "S_ALSFRS_Others", "N_ALSFRS_ALS", "N_ALSFRS_Others", "F_ALSFRS_ALS", "F_ALSFRS_Others"), col = c("black","darkgrey","black","darkgrey","black","darkgrey"), lty = c("solid","solid","dashed","dashed","dotdash","dotdash"), cex = 0.6)
}



