# Simulating data from an AFT model with time-varying covariates and partly-interval censoring:

dataGenAFT <- function(n, beta, gamma, tau_min, tau_max, dist, alpha, psi, pi_E, alpha_L, alpha_R){
  
  if(missing(tau_min)){tau_min = 0}
  if(missing(tau_max)){tau_max = 1}
  if(missing(alpha_L)){alpha_L = 0.5}
  if(missing(alpha_R)){alpha_R = 1.5}
  
  ####################################################
  #### Generate matrix for time-fixed covariates: ####
  ####################################################
  
  x1 = rbinom(n, 1, 0.5)
  x2 = runif(n, 0, 3)
  
  Xmat = cbind(x1, x2)
  
  #####################################
  #### Generate exact event times: ####
  #####################################
  
  # Generate vector of standard uniform random numbers (for inverse transform sampling):
  U = runif(n, 0, 1)
  
  # Generate treatment times:
  tau = runif(n, min = tau_min, max = tau_max)
  
  # Select hazard type for T_0: dist = "weibull" or dist = "log-logistic"
  # Set values for the scale and shape parameter: E.g. alpha = 3 and psi = 1
  
  if(dist=="weibull"){
    # Generate t values according to t < tau:
    t = exp(Xmat%*%beta)*((-(1/psi)*log(U))^(1/alpha)) 
    
    # Check which values satisfy t < tau:
    treat_index = t < tau
    
    # Generate remaining t according to S^{-1}(u) for t > tau:
    t_alt = exp(gamma)*(exp(Xmat%*%beta)*((-(1/psi)*log(U))^(1/alpha)) - tau) + tau
    t[treat_index == FALSE] = t_alt[treat_index == FALSE]
  }
  
  if(dist=="log-logistic"){
    # If t < tau:
    t = exp(Xmat%*%beta)*(((1/psi)*(1/U - 1))^(1/alpha))
    
    # Check which values satisfy t < tau:
    treat_index = t < tau
    
    # Generate remaining t according to S^{-1}(u) for t > tau:
    t_alt = exp(gamma)*(exp(Xmat%*%beta)*(((1/psi)*(1/U - 1))^(1/alpha)) - tau) + tau
    t[treat_index == FALSE] = t_alt[treat_index == FALSE]
  }
  
  ####################################################
  #### Generate vectors for censoring indicators: ####
  ####################################################
  
  # Set proportion of exact event times using
  # E.g. pi_E = 0.3
  
  # Generate standard uniform r.v. for event times:
  UE = runif(n, 0, 1)
  
  # Generate standard uniform r.v. for interval censoring times:
  UL = runif(n, 0, 1)
  UR = runif(n, UL, 1)
  Umat = cbind(UL, UR)
  Umat = apply(Umat, 1, sort)
  UL = Umat[1,]
  UR = Umat[2,]
  
  # Select scalars to define width of censoring intervals:
  
  # Decreasing this lowers proportion of left-censoring & increases proportion of interval censoring:
  # E.g. alpha_L = 0.35
  
  # Increasing this lowers proportion of right-censoring & increases proportion of interval censoring:
  # E.g. alpha_R = 0.75
  
  # Assess range of intervals:
  #range(alpha_R*UR - alpha_L *UL)
  #boxplot(alpha_R*UR - alpha_L *UL)
  
  # Generate the various censoring indicators:
  delE = as.vector((UE < pi_E)*1)
  delL = as.vector(((UE >= pi_E) & (t < alpha_L*UL))*1)
  delR = as.vector(((UE >= pi_E) & (alpha_R*UR < t))*1)
  delI = as.vector(((UE >= pi_E) & (alpha_L*UL <= t) &  (t <= alpha_R*UR))*1)
  
  del = list(delE, delL, delR, delI)
  
  # Check that proportions sum to 1:
  censor_prop = apply(cbind(delE, delL, delR, delI), 2, sum)/n
  names(censor_prop) <- c("Event times", "Left-censored", "Right-censored", "Interval-censored")
  
  #######################################################
  #### Generate vectors for observed survival times: ####
  #######################################################
  yL = (t^(delE))*((alpha_L*UL)^delI)*((alpha_R*UR)^delR)*(0^delL)
  yR = (t^(delE))*((alpha_L*UL)^delL)*((alpha_R*UR)^delI)*(Inf^delR)
  
  # Check that all entries in y^R are >= y^L:
  #sum(yR>=yL) == n
  
  ######################################################
  #### Generate matrix for time-varying covariates: ####
  ######################################################
  # Value of t_ni: max(yL, yR) < Inf:
  t_ni = (yL^delE)*(yR^delL)*(yL^delR)*(yR^delI)
  
  ## Generate time-varying covariates in long-format: ##
  # Create Z matrix for time-varying covariates with a single change-point:
  
  Zmat = {}
  
  for(i in 1:n){
    if(tau[i] > t_ni[i]){
      Zmat = rbind(Zmat, c(i, 0, t_ni[i], 0))
    }
    
    else{
      Zmat = rbind(Zmat, c(i, 0, tau[i], 0), c(i, tau[i], t_ni[i], 1))
    }
  }
  
  colnames(Zmat) = c("Individual", "Start", "End", "z_i(t)")
  
  ## Generate time-varying covariates for z_i(t_ni): ##
  zfin = (t_ni >= tau)*1 
  
  ##################################################
  #### Summarise and return all output objects: ####
  ##################################################
  return(list(exact = t, Xmat = Xmat, Zmat = Zmat, del = del, censor_prop = censor_prop, 
              UL = UL, UR = UR, yL = yL, yR = yR, t_ni = t_ni, zfin = zfin, treat_index = treat_index))
}
