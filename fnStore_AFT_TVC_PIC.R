# Necessary functions for both optimisation and numerical differentiation: 

#########################################################
#### Function to calculate penalized log-likelihood: ####
#########################################################

penloglikfun <- function(beta, gamma, theta, Xmat, knots, sd, h_Rmat,
                         id_Z, Z_t, time_int, mod_time_int,
                         delE, delL, delR, delI, zfin){
  
  # Calculate current kappa values:
  eXBeta = exp(-Xmat%*%beta)
  eZGamma = exp(-Z_t%*%gamma)
    
  kp_right = eXBeta*as.matrix(tapply(eZGamma*time_int, id_Z, sum))
  kp_left_int = eXBeta*as.matrix(tapply(eZGamma*mod_time_int, id_Z, sum))
  
  # Baseline hazard and cumulative baseline hazard (right points):
  kpstd_right = apply(as.matrix(kp_right), 1, function(x){(x-knots)/sd})
  
  cumgsb_right = pnorm(kpstd_right) - pnorm(-knots/sd)
  c_haz = as.vector(theta%*%cumgsb_right)
  
  gsb_right = dnorm(kpstd_right)/sd
  haz = as.vector(theta%*%gsb_right) 
  
  # Cumulative baseline hazard (for left points from interval censoring):
  kpstd_left = apply(as.matrix(kp_left_int), 1, function(x){(x-knots)/sd})
  
  cumgsb_left = pnorm(kpstd_left) - pnorm(-knots/sd)
  c_haz_left = as.vector(theta%*%cumgsb_left)
  
  # Compute the log-likelihood:
  loglikE = delE*(log(haz) - (Xmat%*%beta + zfin*gamma) - c_haz)
  loglikR = delR*c_haz
  loglikL = delL*log1p(-exp(-c_haz))  #loglikL = delL*log(1 - exp(-c_haz))
  loglikI = delI*log((expm1(-c_haz_left) - expm1(-c_haz))^delI) 
  
  loglik = sum(loglikE - loglikR + loglikL + loglikI)
  
  # Compute the penalized log-likelihood:
  penloglik = loglik - theta%*%h_Rmat%*%theta
  
  return(penloglik)
}

####################################################
#### Calculate kappa (accelerated time) values: ####
####################################################

kpvec_fun <- function(beta, gamma, XmatIn, Z_t, id_Z, time_int, mod_time_int){
  
  eXBeta = exp(-XmatIn%*%beta)
  eZGamma = exp(-Z_t%*%gamma)
  
  eZGammaTime = eZGamma*time_int
  eZGammaModTime = eZGamma*mod_time_int
    
  kp_right =  eXBeta*as.matrix(tapply(eZGammaTime, id_Z, sum))
  kp_d_right = -eXBeta*as.matrix(tapply(eZGammaTime*Z_t, id_Z, sum))
  kp_2d_right =  eXBeta*as.matrix(tapply(eZGammaTime*Z_t*Z_t, id_Z, sum))
  
  # Final k_i(y_i^L), for observations with interval censoring:
  kp_left_int = eXBeta*as.matrix(tapply(eZGammaModTime, id_Z, sum))
  kp_d_left_int = -eXBeta*as.matrix(tapply(eZGammaModTime*Z_t, id_Z, sum))
  kp_2d_left_int = eXBeta*as.matrix(tapply(eZGammaModTime*Z_t*Z_t, id_Z, sum))
  
  return(list(kp_right = kp_right, kp_d_right = kp_d_right, kp_2d_right = kp_2d_right,
              kp_left_int = kp_left_int, kp_d_left_int = kp_d_left_int, kp_2d_left_int = kp_2d_left_int))
  
}

#################################
#### Knots and sigma values: ####
#################################
# This function generates the knots and sigma values for the Gaussian basis functions:
# Note: kappa_vec is a vector of k_i(y_i), i = 1,...,n.

knots_fun <- function(num_basis, knots_option, quantSelect, sd_option, kappa_vec){
  
  min_quantile = quantSelect[1]
  max_quantile = quantSelect[2]
  
  if(knots_option == "equal-space"){
    
    # Calculate positions of equally-spaced knots:
    range_kp = range(kappa_vec)
    gknots = seq(range_kp[1], range_kp[2], length.out = num_basis)
    
    # Calculate sigma:
    bin_width = diff(gknots)[1]
    sig = rep(((1/2)*bin_width), num_basis)
    
    return(list(gknots=gknots, sd = sig))
  }
  
  if(knots_option == "quantile"){
    
    # Calculate positions of knots based on quantile approach:
    gknots = quantile(kappa_vec, probs = seq(min_quantile, max_quantile,
                                             length.out = num_basis), names = FALSE)
    
    # Calculate sigma:
    bin_width = diff(gknots)
    
    
    if(sd_option == 1){
      bin_width = c(0.98*bin_width[1], bin_width[1:(num_basis-2)], 0.98*bin_width[num_basis-1])
      sig = c((1/2)*bin_width) # Half of bin width:
    }
    else if(sd_option == 2){
      dist = sapply(gknots, function(x){abs(x - kappa_vec)})
      sig = apply(dist, 2, function(x){quantile(x, 2/3)})/2 # Covers 2/3 of datapoints:
    }
    else{
      sig = NA
      print("sd_option is required.")}
    
    return(list(gknots = gknots, sd = sig))
  }
}

#########################################
#### Hazards and cumulative hazards: ####
#########################################
# This function generates vectors for the hazards and cumulative hazards:

gauss_basis_fun <- function(theta, kappa, knots, sd){
  
  # Returns dimension (number of basis functions (x) n)
  kpstd = apply(as.matrix(kappa), 1, function(x){(x-knots)/sd})
  
  cumgsb = pnorm(kpstd) - pnorm(-knots/sd)
  c_haz_gen = as.vector(theta%*%cumgsb)
  
  gsb = dnorm(kpstd)/sd
  haz_gen = as.vector(theta%*%gsb)
  
  gsb_d = (-1/(sd^2))*kpstd*dnorm(kpstd)
  haz_d_gen = as.vector(theta%*%gsb_d)
  
  gsb_2d = (1/(sd^3))*((kpstd^2)-1)*dnorm(kpstd)
  haz_2d_gen = as.vector(theta%*%gsb_2d)
  
  return(list(c_haz_gen=c_haz_gen, cumgsb = cumgsb,
              haz_gen=haz_gen, haz_d_gen = haz_d_gen, haz_2d_gen = haz_2d_gen, 
              gsb = gsb, gsb_d = gsb_d, gsb_2d = gsb_2d))
}

########################################################
#### Computing gradients for beta, gamma and theta: ####
########################################################
grad_fun <- function(beta, gamma, theta, knots, sd, h_Rmat, XmatIn,
                     id_Z, Z_t, time_int, mod_time_int, delE, delL, delR, delI, zfin){
  
  ### Calculate kappa values: ###
  kpset = kpvec_fun(beta, gamma, XmatIn, Z_t, id_Z, time_int, mod_time_int)
  
  kp_right = kpset$kp_right
  kp_d_right = kpset$kp_d_right
  
  kp_left_int = kpset$kp_left_int
  kp_d_left_int = kpset$kp_d_left_int
  
  ### Calculate hazard values at this iteration: ###
  # Baseline hazard and cumulative baseline hazard (most cases):
  set1 = gauss_basis_fun(theta, kp_right, knots, sd)
  c_haz = set1$c_haz_gen
  haz = set1$haz_gen
  haz_d = set1$haz_d_gen
  gsb_all = t(set1$gsb)
  cumgsb_all = t(set1$cumgsb)
  
  # Cumulative baseline hazard (When considering the left point in int. censoring):
  set2 = gauss_basis_fun(theta, kp_left_int, knots, sd)
  c_haz_left = set2$c_haz_gen
  haz_left = set2$haz_gen
  gsb_left = t(set2$gsb)
  cumgsb_left = t(set2$cumgsb)
  
  # Compute common quantities:
  eCHaz =  exp(-c_haz)
  eCHazLeft = exp(-c_haz_left)
  
  ### Fix numerical issues (for indices associated with event times): ###
  exp2denfix = (eCHazLeft - eCHaz)^delI
  
  ### Calculate gradient of beta at this iteration: ###
  expLall = eCHazLeft*haz_left*kp_left_int
  expRall = eCHaz*haz*kp_right
  
  E_part_gb = delE*((-(haz_d*kp_right)/haz) -1 + haz*kp_right)
  R_part_gb = delR*(haz*kp_right)
  L_part_gb = delL*(expRall/(1-eCHaz))
  I_part_gb = delI*((expLall - expRall)/exp2denfix)
  
  grad_beta = colSums(as.vector(E_part_gb + R_part_gb - L_part_gb + I_part_gb)*XmatIn)
  
  ### Calculate gradient of gamma at this iteration: ###
  expdLall = eCHazLeft*haz_left*kp_d_left_int
  expdRall = eCHaz*haz*kp_d_right
  
  E_part_gg = delE*(-zfin + (haz_d/haz - haz)*kp_d_right)
  R_part_gg = delR*(haz*kp_d_right)
  L_part_gg = delL*((eCHaz*haz)/(1-eCHaz))*kp_d_right
  I_part_gg = delI*(-expdLall + expdRall)/exp2denfix
  
  grad_gamma = colSums(E_part_gg - R_part_gg + L_part_gg + I_part_gg)
  
  ### Calculate gradient of theta at this iteration: ###
  E_part_gt = delE*((gsb_all/haz - cumgsb_all)^delE)
  R_part_gt = delR*((cumgsb_all)^delR)
  L_part_gt = delL*(((cumgsb_all*eCHaz)/(1-eCHaz))^delL)
  I_part_gt = delI*(((-cumgsb_left*eCHazLeft + cumgsb_all*eCHaz)/exp2denfix)^delI)
  
  penalty_d_gt = 2*h_Rmat%*%theta
  
  grad_theta = colSums(E_part_gt - R_part_gt + L_part_gt + I_part_gt) - penalty_d_gt

  return(c(grad_beta, grad_gamma, grad_theta))
}