plAFT_GIC_TVC_Gaussian_1_autosmooth_WBRT_plot <- function(data, estimation_result, tau = 2, resolution = 500, plot_sig_lv = 0.05) {
  
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
  
  gknots_plot = estimation_result$knots_location
  sig_plot = estimation_result$knots_sigma
  cov_mat_theta = estimation_result$cov_mat_theta
  index_active_theta = estimation_result$index_active_theta
  
  Gaussian_basis <- function(x, mean, sd) {
    if (nrow(x)==1) {return(matrix(sapply(X = 1:length(mean), FUN = function(i) {dnorm(x, mean[i], sd[i])}), nrow = 1))}
    else {return(sapply(X = 1:length(mean), FUN = function(i) {dnorm(x, mean[i], sd[i])}))}
  }
  Gaussian_basis_integ1 <- function(x, mean, sd) {
    if (nrow(x)==1) {return(matrix(sapply(X = 1:length(mean), FUN = function(i) {(pnorm(x, mean[i], sd[i]) - matrix(1, NROW(x), 1)%*%pnorm(0, mean[i], sd[i]))}), nrow = 1))}
    else {return(sapply(X = 1:length(mean), FUN = function(i) {(pnorm(x, mean[i], sd[i]) - matrix(1, NROW(x), 1)%*%pnorm(0, mean[i], sd[i]))}))}
  }
  
  t_plot = seq(min(1e-10, range_time[1]), range_time[2], length.out = resolution)
  
  # create the Zmat based on the records of {0, 0}
  Zmat_00 = NULL
  for (i in 1:length(t_plot)) {
    if (t_plot[i] < tau) {
      Zmat_00 = rbind(Zmat_00,
                      c(i, 1, 0, t_plot[i], 0))
    }
    else if (t_plot[i] < tau+3/12) {
      if (tau == 0) {
        Zmat_00 = rbind(Zmat_00,
                        c(i, 1, 0, t_plot[i], 0))
      }
      else {
        Zmat_00 = rbind(Zmat_00,
                        c(i, 1, 0, tau, 0),
                        c(i, 2, tau, t_plot[i], 0))
      }
    }
    else {
      if (tau == 0) {
        Zmat_00 = rbind(Zmat_00,
                        c(i, 1, 0, 3/12, 0),
                        c(i, 2, 3/12, t_plot[i], 0))
      }
      else {
        Zmat_00 = rbind(Zmat_00,
                        c(i, 1, 0, tau, 0),
                        c(i, 2, tau, tau+3/12, 0),
                        c(i, 3, tau+3/12, t_plot[i], 0))
      }
    }
  }
  rownames(Zmat_00) = NULL
  id_Z_00 = Zmat_00[,1]
  interval_length_00 = Zmat_00[, 4, drop = FALSE] - Zmat_00[, 3, drop = FALSE] # N by 1 matrix
  Z_00 = Zmat_00[, -c(1:4), drop = FALSE] # N by q matrix
  Z_00_ni = Z_00[c(diff(id_Z_00)!=0,TRUE), , drop = FALSE]
  
  # create the Zmat based on the records of {0, 1}
  Zmat_01 = NULL
  for (i in 1:length(t_plot)) {
    if (t_plot[i] < tau) {
      Zmat_01 = rbind(Zmat_01,
                      c(i, 1, 0, t_plot[i], 0))
    }
    else if (t_plot[i] < tau+3/12) {
      if (tau == 0) {
        Zmat_01 = rbind(Zmat_01,
                        c(i, 1, 0, t_plot[i], 1))
      }
      else {
        Zmat_01 = rbind(Zmat_01,
                        c(i, 1, 0, tau, 0),
                        c(i, 2, tau, t_plot[i], 1))
      }
    }
    else {
      if (tau == 0) {
        Zmat_01 = rbind(Zmat_01,
                        c(i, 1, 0, 3/12, 1),
                        c(i, 2, 3/12, t_plot[i], 0))
      }
      else {
        Zmat_01 = rbind(Zmat_01,
                        c(i, 1, 0, tau, 0),
                        c(i, 2, tau, tau+3/12, 1),
                        c(i, 3, tau+3/12, t_plot[i], 0))
      }
    }
  }
  rownames(Zmat_01) = NULL
  id_Z_01 = Zmat_01[,1]
  interval_length_01 = Zmat_01[, 4, drop = FALSE] - Zmat_01[, 3, drop = FALSE] # N by 1 matrix
  Z_01 = Zmat_01[, -c(1:4), drop = FALSE] # N by q matrix
  Z_01_ni = Z_01[c(diff(id_Z_01)!=0,TRUE), , drop = FALSE]
  
  
  plot_ingredients <- function(X_input, beta = beta_hat, Z_input, Z_ni_input, gamma = gamma_hat, theta = theta_hat, interval_length_input, id_Z_input, cov_mat_theta_input = cov_mat_theta, index_active_theta_input = index_active_theta) {
    kappa_plot = exp(-X_input%*%beta) * as.matrix(tapply(exp(-Z_input%*%gamma)*interval_length_input, id_Z_input, FUN = sum))
    
    knots = gknots_plot
    sigs = sig_plot
    
    phi_plot = Gaussian_basis(x = kappa_plot, mean = knots, sd = sigs)
    Phi_plot = Gaussian_basis_integ1(x = kappa_plot, mean = knots, sd = sigs)
    
    if (length(index_active_theta_input)>0) {
      cov_mat_theta_no_active = cov_mat_theta_input[-index_active_theta_input,-index_active_theta_input]
      phi_plot_no_active = phi_plot[,-index_active_theta_input]
      h0kappa_plot_no_active = phi_plot_no_active%*%theta[-index_active_theta_input]
      Phi_plot_no_active = Phi_plot[,-index_active_theta_input]
      S0kappa_plot_no_active = exp(-Phi_plot_no_active%*%theta[-index_active_theta_input])
    }
    else {
      cov_mat_theta_no_active = cov_mat_theta_input
      phi_plot_no_active = phi_plot
      h0kappa_plot_no_active = phi_plot%*%theta
      Phi_plot_no_active = Phi_plot
      S0kappa_plot_no_active = exp(-Phi_plot%*%theta)
    }
    
    # baseline hazard plot of accelerated/decelerated time scale
    G_h0kappa_plot = log(h0kappa_plot_no_active)
    dG_dtheta_h0kappa_plot = (1/h0kappa_plot_no_active)%*%matrix(1,1,ncol(phi_plot_no_active)) * phi_plot_no_active
    asym_std_G_h0kappa_plot = sqrt(diag(dG_dtheta_h0kappa_plot%*%cov_mat_theta_no_active%*%t(dG_dtheta_h0kappa_plot)))
    CI_G_h0kappa_plot = cbind(G_h0kappa_plot-qnorm(1-plot_sig_lv/2)*asym_std_G_h0kappa_plot, G_h0kappa_plot+qnorm(1-plot_sig_lv/2)*asym_std_G_h0kappa_plot)
    CI_h0kappa_plot = exp(CI_G_h0kappa_plot)
    
    ht_plot = h0kappa_plot_no_active * exp(-X_input%*%beta - Z_ni_input%*%gamma)
    G_ht_plot = log(ht_plot)
    asym_std_G_ht_plot = asym_std_G_h0kappa_plot
    CI_G_ht_plot = cbind(G_ht_plot-qnorm(1-plot_sig_lv/2)*asym_std_G_ht_plot, G_ht_plot+qnorm(1-plot_sig_lv/2)*asym_std_G_ht_plot)
    CI_ht_plot = exp(CI_G_ht_plot)
    
    # baseline survival plot of accelerated/decelerated time scale
    G_S0kappa_plot = log(S0kappa_plot_no_active/(1-S0kappa_plot_no_active))
    dG_dtheta_S0kappa_plot = -(1/(1-S0kappa_plot_no_active))%*%matrix(1,1,ncol(Phi_plot_no_active)) * Phi_plot_no_active
    asym_std_G_S0kappa_plot = sqrt(diag(dG_dtheta_S0kappa_plot%*%cov_mat_theta_no_active%*%t(dG_dtheta_S0kappa_plot)))
    CI_G_S0kappa_plot = cbind(G_S0kappa_plot-qnorm(1-plot_sig_lv/2)*asym_std_G_S0kappa_plot, G_S0kappa_plot+qnorm(1-plot_sig_lv/2)*asym_std_G_S0kappa_plot)
    CI_S0kappa_plot = 1/(exp(-CI_G_S0kappa_plot)+1)
    
    results_plot_ingredients = list()
    results_plot_ingredients$X_input = X_input
    results_plot_ingredients$kappa_plot = kappa_plot
    
    results_plot_ingredients$theta_hat = theta_hat
    results_plot_ingredients$phi_plot = phi_plot
    
    results_plot_ingredients$h0kappa_plot = h0kappa_plot_no_active
    results_plot_ingredients$CI_h0kappa_plot = CI_h0kappa_plot
    
    results_plot_ingredients$ht_plot = ht_plot
    results_plot_ingredients$CI_ht_plot = CI_ht_plot
    
    results_plot_ingredients$S0kappa_plot = S0kappa_plot_no_active
    results_plot_ingredients$CI_S0kappa_plot = CI_S0kappa_plot
    
    return(results_plot_ingredients)
  }
  
  # design X input
  # Plot for Treat: Treat = 1, Gender = 0, Number = 1, Extracranial_at_Baseline = 1, Age = median(Age)
  X_T1 = matrix(1, resolution, 5)
  X_T1[, 2] = 0
  X_T1[, 5] = median(X[, 5])
  
  option_X_T1_00 = plot_ingredients(X_input = X_T1, Z_input = Z_00, Z_ni_input = Z_00_ni, interval_length_input = interval_length_00, id_Z_input = id_Z_00)
  option_X_T1_01 = plot_ingredients(X_input = X_T1, Z_input = Z_01, Z_ni_input = Z_01_ni, interval_length_input = interval_length_01, id_Z_input = id_Z_01)
  
  option_x_T1_00_tau = plot_ingredients(X_input = X_T1[1:2, , drop = FALSE], Z_input = Z_00[1:2, , drop = FALSE], Z_ni_input = Z_00_ni[1:2, , drop = FALSE], interval_length_input = matrix(tau, 2, 1), id_Z_input = c(1,2))
  option_x_T1_01_tau = plot_ingredients(X_input = X_T1[1:2, , drop = FALSE], Z_input = Z_01[1:2, , drop = FALSE], Z_ni_input = Z_01_ni[1:2, , drop = FALSE], interval_length_input = matrix(tau, 2, 1), id_Z_input = c(1,2))
  
  # Plot for Treat: Treat = 0, Gender = 0, Number = 1, Extracranial_at_Baseline = 1, Age = median(Age)
  X_T0 = matrix(1, resolution, 5)
  X_T0[, 1:2] = 0
  X_T0[, 5] = median(X[, 5])
  
  option_X_T0_00 = plot_ingredients(X_input = X_T0, Z_input = Z_00, Z_ni_input = Z_00_ni, interval_length_input = interval_length_00, id_Z_input = id_Z_00)
  option_X_T0_01 = plot_ingredients(X_input = X_T0, Z_input = Z_01, Z_ni_input = Z_01_ni, interval_length_input = interval_length_01, id_Z_input = id_Z_01)
  
  option_x_T0_00_tau = plot_ingredients(X_input = X_T0[1:2, , drop = FALSE], Z_input = Z_00[1:2, , drop = FALSE], Z_ni_input = Z_00_ni[1:2, , drop = FALSE], interval_length_input = matrix(tau, 2, 1), id_Z_input = c(1,2))
  option_x_T0_01_tau = plot_ingredients(X_input = X_T0[1:2, , drop = FALSE], Z_input = Z_01[1:2, , drop = FALSE], Z_ni_input = Z_01_ni[1:2, , drop = FALSE], interval_length_input = matrix(tau, 2, 1), id_Z_input = c(1,2))
  
  
  # assuming patient has survived up to tau years, what will be the probability of survival later on
  index_tau = min(which(t_plot > tau))
  dp_time_plot_tau = c(0, tau, t_plot[index_tau:resolution])
  # X_T1_00 and X_T1_01
  dp_option_X_T1_00_tau = c(1, 1, option_X_T1_00$S0kappa_plot[index_tau:resolution]/option_x_T1_00_tau$S0kappa_plot[1])
  dp_option_X_T1_01_tau = c(1, 1, option_X_T1_01$S0kappa_plot[index_tau:resolution]/option_x_T1_01_tau$S0kappa_plot[1])
  # X_T0_00 and X_T0_01
  dp_option_X_T0_00_tau = c(1, 1, option_X_T0_00$S0kappa_plot[index_tau:resolution]/option_x_T0_00_tau$S0kappa_plot[1])
  dp_option_X_T0_01_tau = c(1, 1, option_X_T0_01$S0kappa_plot[index_tau:resolution]/option_x_T0_01_tau$S0kappa_plot[1])
  
  # plot
  plot(dp_time_plot_tau, dp_option_X_T0_00_tau, xlim = c(0,2.5), ylim = c(0.3,1), xlab = expression(t), ylab = bquote(P(T[i]~">"~tau[i]+dt~"|"~T[i]~">"~tau[i],bold(x)[i],tilde(bold(z))[i](t),tau[i]==.(tau))), mgp=c(2,1,0), type = "l", lty = "solid", col = "white") # 20241030
  grid(lty = "solid")
  lines(dp_time_plot_tau, dp_option_X_T1_00_tau, type = "l", lty = "dotdash", lwd = 2, col = "black")
  lines(dp_time_plot_tau, dp_option_X_T1_01_tau, type = "l", lty = "solid", lwd = 2, col = "black")
  lines(dp_time_plot_tau, dp_option_X_T0_00_tau, type = "l", lty = "dotdash", lwd = 2, col = "blue")
  lines(dp_time_plot_tau, dp_option_X_T0_01_tau, type = "l", lty = "solid", lwd = 2, col = "blue")
  abline(v = tau, lty = "dashed", lwd = 2, col = "red")
  text(tau+0.1, 0.65, bquote(tau[i]==.(tau)), srt = 90, col = "red", cex = 0.8)
  title(main = paste("survived beyond", tau, "year(s)"), cex = 0.8)
  legend("bottomright", legend = c("With WBRT, systemic therapy: yes", "With WBRT, systemic therapy: no", "Without WBRT, systemic therapy: yes", "Without WBRT, systemic therapy: no"), text.width = strwidth("Without WBRT, systemic therapy: yes")[1], col = c("black","black","blue","blue"), lty = c("solid","dotdash","solid","dotdash"), lwd = 2, cex = 0.8)
}



