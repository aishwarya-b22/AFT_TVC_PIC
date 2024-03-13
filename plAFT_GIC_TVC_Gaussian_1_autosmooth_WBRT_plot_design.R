# forked from plAFT_GIC_TVC_Gaussian_1_autosmooth_MND_plot_design.R on 20240118
plAFT_GIC_TVC_Gaussian_1_autosmooth_WBRT_plot_design <- function(data, estimation_result, covariate = "Treat", resolution = 500, knots_option = "percentile", plot_sig_lv = 0.05) {
  # browser()
  tmat = as.matrix(data$tmat)
  rownames(tmat) = NULL
  time = c(tmat[,1],tmat[,2])
  time = time[!is.infinite(time)]
  range_time = range(time)
  censor_n = tmat[, 3, drop = FALSE] # n by 1 matrix
  
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
    if (time_plot[i] < 2) {
      Zmat_00 = rbind(Zmat_00,
                      c(i, 1, 0, time_plot[i], 0))
    }
    else {
      Zmat_00 = rbind(Zmat_00,
                      c(i, 1, 0, 2, 0),
                      c(i, 2, 2, time_plot[i], 0))
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
    if (time_plot[i] < 2) {
      Zmat_01 = rbind(Zmat_01,
                      c(i, 1, 0, time_plot[i], 0))
    }
    else {
      Zmat_01 = rbind(Zmat_01,
                      c(i, 1, 0, 2, 0),
                      c(i, 2, 2, time_plot[i], 1))
    }
  }
  rownames(Zmat_01) = NULL
  id_Z_01 = Zmat_01[,1]
  interval_length_01 = Zmat_01[, 4, drop = FALSE] - Zmat_01[, 3, drop = FALSE] # N by 1 matrix
  Z_01 = Zmat_01[, -c(1:4), drop = FALSE] # N by q matrix
  Z_01_ni = Z_01[c(diff(id_Z_01)!=0,TRUE), , drop = FALSE]
  
  # create the Zmat based on the records of {1, 0}
  Zmat_10 = NULL
  for (i in 1:length(time_plot)) {
    if (time_plot[i] < 2) {
      Zmat_10 = rbind(Zmat_10,
                      c(i, 1, 0, time_plot[i], 1))
    }
    else {
      Zmat_10 = rbind(Zmat_10,
                      c(i, 1, 0, 2, 1),
                      c(i, 2, 2, time_plot[i], 0))
    }
  }
  rownames(Zmat_10) = NULL
  id_Z_10 = Zmat_10[,1]
  interval_length_10 = Zmat_10[, 4, drop = FALSE] - Zmat_10[, 3, drop = FALSE] # N by 1 matrix
  Z_10 = Zmat_10[, -c(1:4), drop = FALSE] # N by q matrix
  Z_10_ni = Z_10[c(diff(id_Z_10)!=0,TRUE), , drop = FALSE]
  
  # create the Zmat based on the records of {1, 1}
  Zmat_11 = NULL
  for (i in 1:length(time_plot)) {
    if (time_plot[i] < 2) {
      Zmat_11 = rbind(Zmat_11,
                      c(i, 1, 0, time_plot[i], 1))
    }
    else {
      Zmat_11 = rbind(Zmat_11,
                      c(i, 1, 0, 2, 1),
                      c(i, 2, 2, time_plot[i], 1))
    }
  }
  rownames(Zmat_11) = NULL
  id_Z_11 = Zmat_11[,1]
  interval_length_11 = Zmat_11[, 4, drop = FALSE] - Zmat_11[, 3, drop = FALSE] # N by 1 matrix
  Z_11 = Zmat_11[, -c(1:4), drop = FALSE] # N by q matrix
  Z_11_ni = Z_11[c(diff(id_Z_11)!=0,TRUE), , drop = FALSE]
  
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
    X_T1 = matrix(1, resolution, 5)
    if (covariate == "Treat") {# Plot for Treat: Treat = 1, Gender = 0, Number = 1, Extracranial_at_Baseline = 1, Age = median(Age)
      X_T1[, 2] = 0
      X_T1[, 5] = median(X[, 5])
    }
    if (covariate == "Gender") {# Plot for Gender: Treat = 1, Gender = 1, Number = 1, Extracranial_at_Baseline = 1, Age = median(Age)
      X_T1[, 5] = median(X[, 5])
    }
    if (covariate == "Number") {# Plot for Number: Treat = 1, Gender = 0, Number = 1, Extracranial_at_Baseline = 1, Age = median(Age)
      X_T1[, 2] = 0
      X_T1[, 5] = median(X[, 5])
    }
    if (covariate == "ED") {# Plot for Extracranial_at_Baseline: Treat = 1, Gender = 0, Number = 1, Extracranial_at_Baseline = 1, Age = median(Age)
      X_T1[, 2] = 0
      X_T1[, 5] = median(X[, 5])
    }
    if (covariate == "Age") {# Plot for Age groups (high): Treat = 1, Gender = 0, Number = 1, Extracranial_at_Baseline = 1, Age = median(higher age group)
      X_T1[, 2] = 0
      X_T1[, 5] = median(X[X[, 5] >= quantile(X[, 5], seq(0, 1, length.out = 4))[3], 5])
    }
    option_X_T1_00 = plot_ingredients(X_input = X_T1, Z_input = Z_00, Z_ni_input = Z_00_ni, interval_length_input = interval_length_00, id_Z_input = id_Z_00)
    option_X_T1_01 = plot_ingredients(X_input = X_T1, Z_input = Z_01, Z_ni_input = Z_01_ni, interval_length_input = interval_length_01, id_Z_input = id_Z_01)
    option_X_T1_10 = plot_ingredients(X_input = X_T1, Z_input = Z_10, Z_ni_input = Z_10_ni, interval_length_input = interval_length_10, id_Z_input = id_Z_10)
    option_X_T1_11 = plot_ingredients(X_input = X_T1, Z_input = Z_11, Z_ni_input = Z_11_ni, interval_length_input = interval_length_11, id_Z_input = id_Z_11)
  }
  
  {
    X_T0 = matrix(1, resolution, 5)
    if (covariate == "Treat") {# Plot for Treat: Treat = 0, Gender = 0, Number = 1, Extracranial_at_Baseline = 1, Age = median(Age)
      X_T0[, 1:2] = 0
      X_T0[, 5] = median(X[, 5])
    }
    if (covariate == "Gender") {# Plot for Gender: Treat = 1, Gender = 0, Number = 1, Extracranial_at_Baseline = 1, Age = median(Age)
      X_T0[, 2] = 0
      X_T0[, 5] = median(X[, 5])
    }
    if (covariate == "Number") {# Plot for Number: Treat = 1, Gender = 0, Number = 0, Extracranial_at_Baseline = 1, Age = median(Age)
      X_T0[, 2:3] = 0
      X_T0[, 5] = median(X[, 5])
    }
    if (covariate == "ED") {# Plot for Extracranial_at_Baseline: Treat = 1, Gender = 0, Number = 1, Extracranial_at_Baseline = 0, Age = median(Age)
      X_T0[, c(2,4)] = 0
      X_T0[, 5] = median(X[, 5])
    }
    if (covariate == "Age") {# Plot for Age groups (middle): Treat = 1, Gender = 0, Number = 1, Extracranial_at_Baseline = 1, Age = median(middle age group)
      X_T0[, 2] = 0
      X_T0[, 5] = median(X[X[, 5] >= quantile(X[, 5], seq(0, 1, length.out = 4))[2] & X[, 5] < quantile(X[, 5], seq(0, 1, length.out = 4))[3], 5])
    }
    option_X_T0_00 = plot_ingredients(X_input = X_T0, Z_input = Z_00, Z_ni_input = Z_00_ni, interval_length_input = interval_length_00, id_Z_input = id_Z_00)
    option_X_T0_01 = plot_ingredients(X_input = X_T0, Z_input = Z_01, Z_ni_input = Z_01_ni, interval_length_input = interval_length_01, id_Z_input = id_Z_01)
    option_X_T0_10 = plot_ingredients(X_input = X_T0, Z_input = Z_10, Z_ni_input = Z_10_ni, interval_length_input = interval_length_10, id_Z_input = id_Z_10)
    option_X_T0_11 = plot_ingredients(X_input = X_T0, Z_input = Z_11, Z_ni_input = Z_11_ni, interval_length_input = interval_length_11, id_Z_input = id_Z_11)
  }
  
  if (covariate == "Age") {# Plot for Age groups (young): Treat = 1, Gender = 0, Number = 1, Extracranial_at_Baseline = 1, Age = median(young age group)
    X_Ty = matrix(1, resolution, 5)
    X_Ty[, 2] = 0
    X_Ty[, 5] = median(X[X[, 5] < quantile(X[, 5], seq(0, 1, length.out = 4))[2], 5])
    option_X_Ty_00 = plot_ingredients(X_input = X_Ty, Z_input = Z_00, Z_ni_input = Z_00_ni, interval_length_input = interval_length_00, id_Z_input = id_Z_00)
    option_X_Ty_01 = plot_ingredients(X_input = X_Ty, Z_input = Z_01, Z_ni_input = Z_01_ni, interval_length_input = interval_length_01, id_Z_input = id_Z_01)
    option_X_Ty_10 = plot_ingredients(X_input = X_Ty, Z_input = Z_10, Z_ni_input = Z_10_ni, interval_length_input = interval_length_10, id_Z_input = id_Z_10)
    option_X_Ty_11 = plot_ingredients(X_input = X_Ty, Z_input = Z_11, Z_ni_input = Z_11_ni, interval_length_input = interval_length_11, id_Z_input = id_Z_11)
  }
  
  # plot(time_plot, option_X_T1_00$S0ts_plot, xlim = c(0,range_time[2]*1.05), ylim = c(0,1), main = "Predicted survival (Treatment)", xlab = expression(t), ylab = expression(S(t)), mgp=c(2,1,0), type = "l", lty = "solid", col = "white") # variables for the paper
  plot(time_plot, option_X_T1_00$S0ts_plot, xlim = c(0,range_time[2]*1.05), ylim = c(0,1), xlab = expression(t), ylab = expression(S(t)), mgp=c(2,1,0), type = "l", lty = "solid", col = "white")
  grid(lty = "solid")
  lines(time_plot, option_X_T1_00$S0ts_plot, type = "l", lty = "solid", lwd = 1.5, col = "black")
  lines(time_plot, option_X_T0_00$S0ts_plot, type = "l", lty = "solid", lwd = 1.5, col = "darkgrey")
  if (covariate == "Age") {lines(time_plot, option_X_Ty_00$S0ts_plot, type = "l", lty = "solid", lwd = 1.5, col = "lightgrey")}
  
  lines(time_plot, option_X_T1_01$S0ts_plot, type = "l", lty = "dashed", lwd = 1.5, col = "black")
  lines(time_plot, option_X_T0_01$S0ts_plot, type = "l", lty = "dashed", lwd = 1.5, col = "darkgrey")
  if (covariate == "Age") {lines(time_plot, option_X_Ty_01$S0ts_plot, type = "l", lty = "dashed", lwd = 1.5, col = "lightgrey")}
  
  lines(time_plot, option_X_T1_10$S0ts_plot, type = "l", lty = "dotdash", lwd = 1.5, col = "black")
  lines(time_plot, option_X_T0_10$S0ts_plot, type = "l", lty = "dotdash", lwd = 1.5, col = "darkgrey")
  if (covariate == "Age") {lines(time_plot, option_X_Ty_10$S0ts_plot, type = "l", lty = "dotdash", lwd = 1.5, col = "lightgrey")}
  
  lines(time_plot, option_X_T1_11$S0ts_plot, type = "l", lty = "dotted", lwd = 1.5, col = "black")
  lines(time_plot, option_X_T0_11$S0ts_plot, type = "l", lty = "dotted", lwd = 1.5, col = "darkgrey")
  if (covariate == "Age") {lines(time_plot, option_X_Ty_11$S0ts_plot, type = "l", lty = "dotted", lwd = 1.5, col = "lightgrey")}
  
  abline(v = 2, lty = "dashed", lwd = 2, col = "red")
  # if (covariate == "Number") {text(3.99-0.2, 0.3, expression(tau==3.99), srt = 90, col = "red", cex = 0.85)}
  # else {text(3.99-0.2, 0.1, expression(tau==3.99), srt = 90, col = "red", cex = 0.85)}
  text(2-0.3, 0.1, expression(tau==2), srt = 90, col = "red", cex = 0.8)
  abline(h = 0.5, lty = "dashed", lwd = 2, col = "blue")
  text(0.75, 0.5-0.3/8, expression(S(t)==0.5), col = "blue", cex = 0.8)
  # legend for Treat
  if (covariate == "Treat") {
    title(main = "Predicted Survival (Treatment)")
    legend("topright", legend = c("00_Treat_WBRT", "00_Treat_Observation", "01_Treat_WBRT", "01_Treat_Observation", "10_Treat_WBRT", "10_Treat_Observation","11_Treat_WBRT", "11_Treat_Observation"), col = c("black","darkgrey","black","darkgrey","black","darkgrey","black","darkgrey"), lty = c("solid","solid","dashed","dashed","dotdash","dotdash","dotted","dotted"), lwd = 2, cex = 0.6)
  }
  if (covariate == "Treat") {
    # 00 plot
    plot(time_plot, option_X_T1_00$S0ts_plot, xlim = c(0,range_time[2]*1.05), ylim = c(0,1), xlab = expression(t), ylab = expression(S(t)), mgp=c(2,1,0), type = "l", lty = "solid", col = "white")
    grid(lty = "solid")
    lines(time_plot, option_X_T1_00$S0ts_plot, type = "l", lty = "solid", lwd = 2, col = "black")
    lines(time_plot, option_X_T0_00$S0ts_plot, type = "l", lty = "dotdash", lwd = 2, col = "black")
    abline(v = 2, lty = "dashed", lwd = 2, col = "red")
    text(2-0.3, 0.1, expression(tau==2), srt = 90, col = "red", cex = 0.8)
    abline(h = 0.5, lty = "dashed", lwd = 2, col = "blue")
    text(0.75, 0.5-0.3/8, expression(S(t)==0.5), col = "blue", cex = 0.8)
    title(main = "00 Predicted Survival (Treatment)")
    legend("topright", legend = c("WBRT", "Observation"), col = c("black","black"), lty = c("solid","dotdash"), lwd = 2, cex = 0.8)
    # 01 plot
    plot(time_plot, option_X_T1_00$S0ts_plot, xlim = c(0,range_time[2]*1.05), ylim = c(0,1), xlab = expression(t), ylab = expression(S(t)), mgp=c(2,1,0), type = "l", lty = "solid", col = "white")
    grid(lty = "solid")
    lines(time_plot, option_X_T1_01$S0ts_plot, type = "l", lty = "solid", lwd = 2, col = "black")
    lines(time_plot, option_X_T0_01$S0ts_plot, type = "l", lty = "dotdash", lwd = 2, col = "black")
    abline(v = 2, lty = "dashed", lwd = 2, col = "red")
    text(2-0.3, 0.1, expression(tau==2), srt = 90, col = "red", cex = 0.8)
    abline(h = 0.5, lty = "dashed", lwd = 2, col = "blue")
    text(0.75, 0.5-0.3/8, expression(S(t)==0.5), col = "blue", cex = 0.8)
    title(main = "01 Predicted Survival (Treatment)")
    legend("topright", legend = c("WBRT", "Observation"), col = c("black","black"), lty = c("solid","dotdash"), lwd = 2, cex = 0.8)
    # 10 plot
    plot(time_plot, option_X_T1_00$S0ts_plot, xlim = c(0,range_time[2]*1.05), ylim = c(0,1), xlab = expression(t), ylab = expression(S(t)), mgp=c(2,1,0), type = "l", lty = "solid", col = "white")
    grid(lty = "solid")
    lines(time_plot, option_X_T1_10$S0ts_plot, type = "l", lty = "solid", lwd = 2, col = "black")
    lines(time_plot, option_X_T0_10$S0ts_plot, type = "l", lty = "dotdash", lwd = 2, col = "black")
    abline(v = 2, lty = "dashed", lwd = 2, col = "red")
    text(2-0.3, 0.1, expression(tau==2), srt = 90, col = "red", cex = 0.8)
    abline(h = 0.5, lty = "dashed", lwd = 2, col = "blue")
    text(0.75, 0.5-0.3/8, expression(S(t)==0.5), col = "blue", cex = 0.8)
    title(main = "10 Predicted Survival (Treatment)")
    legend("topright", legend = c("WBRT", "Observation"), col = c("black","black"), lty = c("solid","dotdash"), lwd = 2, cex = 0.8)
    # 11 plot
    plot(time_plot, option_X_T1_00$S0ts_plot, xlim = c(0,range_time[2]*1.05), ylim = c(0,1), xlab = expression(t), ylab = expression(S(t)), mgp=c(2,1,0), type = "l", lty = "solid", col = "white")
    grid(lty = "solid")
    lines(time_plot, option_X_T1_11$S0ts_plot, type = "l", lty = "solid", lwd = 2, col = "black")
    lines(time_plot, option_X_T0_11$S0ts_plot, type = "l", lty = "dotdash", lwd = 2, col = "black")
    abline(v = 2, lty = "dashed", lwd = 2, col = "red")
    text(2-0.3, 0.1, expression(tau==2), srt = 90, col = "red", lwd = 2, cex = 0.8)
    abline(h = 0.5, lty = "dashed", lwd = 2, col = "blue")
    text(0.75, 0.5-0.3/8, expression(S(t)==0.5), col = "blue", cex = 0.8)
    title(main = "11 Predicted Survival (Treatment)")
    legend("topright", legend = c("WBRT", "Observation"), col = c("black","black"), lty = c("solid","dotdash"), lwd = 2, cex = 0.8)
  }
  
  # legend for Gender
  if (covariate == "Gender") {
    title(main = "Predicted Survival (Gender)")
    legend("topright", legend = c("00_Gender_Male", "00_Gender_Female", "01_Gender_Male", "01_Gender_Female", "10_Gender_Male", "10_Gender_Female","11_Gender_Male", "11_Gender_Female"), col = c("black","darkgrey","black","darkgrey","black","darkgrey","black","darkgrey"), lty = c("solid","solid","dashed","dashed","dotdash","dotdash","dotted","dotted"), lwd = 2, cex = 0.6)
  }
  # legend for Number
  if (covariate == "Number") {
    title(main = "Predicted Survival (Number of MBM)")
    legend("topright", legend = c("00_Num_>1", "00_Num_=1", "01_Num_>1", "01_Num_=1", "10_Num_>1", "10_Num_=1","11_Num_>1", "11_Num_=1"), col = c("black","darkgrey","black","darkgrey","black","darkgrey","black","darkgrey"), lty = c("solid","solid","dashed","dashed","dotdash","dotdash","dotted","dotted"), lwd = 2, cex = 0.6)
  }
  # legend for ED
  if (covariate == "ED") {
    title(main = "Predicted Survival (Extracranial Disease)")
    legend("topright", legend = c("00_ED_Present", "00_ED_Absent", "01_ED_Present", "01_ED_Absent", "10_ED_Present", "10_ED_Absent","11_ED_Present", "11_ED_Absent"), col = c("black","darkgrey","black","darkgrey","black","darkgrey","black","darkgrey"), lty = c("solid","solid","dashed","dashed","dotdash","dotdash","dotted","dotted"), lwd = 2, cex = 0.6)
  }
  # legend for Age
  if (covariate == "Age") {
    title(main = "Predicted Survival (Age Groups)")
    legend("topright", legend = c("00_Age_High", "00_Age_Middle", "00_Age_Low", "01_Age_High", "01_Age_Middle", "01_Age_Low", "10_Age_High", "10_Age_Middle", "10_Age_Low", "11_Age_High", "11_Age_Middle", "11_Age_Low"), col = c("black","darkgrey","lightgrey","black","darkgrey","lightgrey","black","darkgrey","lightgrey","black","darkgrey","lightgrey"), lty = c("solid","solid","solid","dashed","dashed","dashed","dotdash","dotdash","dotdash","dotted","dotted","dotted"), lwd = 2, cex = 0.5)
  }
  
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
  # legend("topright", legend = c("S_ALSFRS_ALS", "S_ALSFRS_Others", "N_ALSFRS_ALS", "N_ALSFRS_Others", "F_ALSFRS_ALS", "F_ALSFRS_Others"), col = c("black","black","blue","blue","red","red"), lty = c("solid","dashed","solid","dashed","solid","dashed"), cex = 0.6)
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
  # legend("topright", legend = c("S_ALSFRS_ALS", "S_ALSFRS_Others", "N_ALSFRS_ALS", "N_ALSFRS_Others", "F_ALSFRS_ALS", "F_ALSFRS_Others"), col = c("black","darkgrey","black","darkgrey","black","darkgrey"), lty = c("solid","solid","dashed","dashed","dotdash","dotdash"), cex = 0.6)
}



