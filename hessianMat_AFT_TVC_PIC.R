# Function to compute the full Hessian for both optimisation: 

fullHessian <- function(beta, gamma, theta, knots, sd, XmatHess, Z_t,
                        id_Z, time_int, mod_time_int, delE, delL, delR, delI){
    
  # Calculate kappa values
  kpset = kpvec_fun(beta, gamma, XmatHess, Z_t, id_Z, time_int, mod_time_int)
   
  kp_right = kpset$kp_right
  kp_d_right = kpset$kp_d_right
  kp_2d_right = kpset$kp_2d_right
  
  kp_left_int = kpset$kp_left_int
  kp_d_left_int = kpset$kp_d_left_int
  kp_2d_left_int = kpset$kp_2d_left_int
  
  # Baseline hazard and cumulative baseline hazard (most cases):
  set1 = gauss_basis_fun(theta, kp_right, knots, sd)
  
  haz = set1$haz_gen
  c_haz = set1$c_haz_gen
  haz_d = set1$haz_d_gen
  haz_2d = set1$haz_2d_gen
  
  gsb_all = set1$gsb
  cumgsb_all = set1$cumgsb
  gsb_d_all = set1$gsb_d
  
  # Cumulative baseline hazard (When considering the left point in int. censoring):
  set2 = gauss_basis_fun(theta, kp_left_int, knots, sd)
  
  haz_left = set2$haz_gen
  c_haz_left = set2$c_haz_gen
  haz_d_left = set2$haz_d_gen
  
  gsb_left = set2$gsb
  cumgsb_left = set2$cumgsb
  
  # Compute common quantities:
  eCHaz =  exp(-c_haz)
  eCHazLeft = exp(-c_haz_left)
  exp1den = 1 - eCHaz 
  exp2denfix = (eCHazLeft - eCHaz)^delI
  
  ######################
  ## Hessian of Beta: ##
  ######################
  
  num1_hb = eCHaz*(haz^2 - haz_d)
  num2_hb = eCHaz*haz
  num3_hb = eCHazLeft*(haz_left^2 - haz_d_left)
  num4_hb = eCHazLeft*haz_left
  
  comp1_hb = num1_hb/exp1den 
  comp2_hb = num2_hb/exp1den
  comp3_hb = num3_hb/exp2denfix
  comp4_hb = num4_hb/exp2denfix
  comp5_hb = num1_hb/exp2denfix
  comp6_hb = num2_hb/exp2denfix
  comp7_hb = (num4_hb*kp_left_int*num2_hb*kp_right)/(exp2denfix^2)
  
  E_part_hb = delE*((((haz*haz_2d)- haz_d^2)/(haz^2) - haz_d)*(kp_right^2) 
                    + ((haz_d/haz) - haz)*kp_right)
  
  R_part_hb = delR*((haz_d*kp_right + haz)*kp_right)
  
  L_part_hb = delL*((comp1_hb + comp2_hb^2)*(kp_right^2) - comp2_hb*kp_right)
  
  I_part_hb = delI*((comp3_hb - comp4_hb^2)*(kp_left_int^2) 
                    - comp4_hb*kp_left_int - (comp5_hb + comp6_hb^2)*(kp_right^2) 
                    + comp6_hb*kp_right + 2*comp7_hb)
  
  hess_beta = t(XmatHess)%*%(as.vector(E_part_hb - R_part_hb - L_part_hb 
                                   + I_part_hb)*XmatHess)
  
  #######################
  ## Hessian of Gamma: ##
  #######################
  num1_hg = eCHaz*haz
  num2_hg = eCHaz*haz^2 - eCHaz*haz_d
  num3_hg = eCHazLeft*haz_left
  num4_hg = eCHazLeft*(haz_left^2) - eCHazLeft*haz_d_left
  
  comp1_hg = num4_hg/exp2denfix
  comp2_hg = num3_hg/exp2denfix
  comp3_hg = num2_hg/exp2denfix
  comp4_hg = num1_hg/exp2denfix
  comp5_hg = (num1_hg*num3_hg)/(exp2denfix^2)
  
  fin1_hg = kp_d_left_int*(comp1_hg - comp2_hg^2)*kp_d_left_int
  fin2_hg = comp2_hg*kp_2d_left_int
  fin3_hg = kp_d_right*(comp3_hg + comp4_hg^2)*kp_d_right
  fin4_hg = comp4_hg*kp_2d_right
  fin5_hg = kp_d_left_int*comp5_hg*kp_d_right
  fin6_hg = kp_d_right*comp5_hg*kp_d_left_int
  
  E_part_hg = delE*(kp_d_right*((haz*haz_2d - haz_d^2)/(haz^2)- haz_d)*kp_d_right 
                 + (haz_d/haz - haz)*kp_2d_right)
  
  R_part_hg = delR*(haz*kp_2d_right + kp_d_right*haz_d*kp_d_right)
  
  L_part_hg = delL*(((eCHaz*haz)/exp1den)*kp_2d_right + 
                   kp_d_right*((eCHaz*haz_d - eCHaz*(haz^2))/exp1den
                             - ((eCHaz*haz)/exp1den)^2)*kp_d_right)
  
  I_part_hg = delI*(fin1_hg - fin2_hg - fin3_hg + fin4_hg + fin5_hg + fin6_hg)
  
  hess_gamma = colSums(E_part_hg - R_part_hg + L_part_hg + I_part_hg)
  
  #######################
  ## Hessian of Theta: ##
  #######################
  
  E_part_ht = -gsb_all%*%(t(gsb_all)*(delE/(haz^2)))
  R_part_ht = 0
  L_part_ht = -cumgsb_all%*%(t(cumgsb_all)*((delL*eCHaz)/((1-eCHaz)^2)))
  
  I_part_1_ht = cumgsb_left%*%(t(cumgsb_left)*((delI*eCHazLeft)/exp2denfix))
  I_part_2_ht = -cumgsb_all%*%(t(cumgsb_all)*((delI*eCHaz)/exp2denfix))  
  I_part_3_ht = -cumgsb_left%*%(t(cumgsb_left)*((delI*(eCHazLeft^2))/(exp2denfix^2)))
  I_part_4_ht = cumgsb_all%*%(t(cumgsb_left)*((delI*(eCHaz*eCHazLeft))/(exp2denfix^2)))
  I_part_5_ht = cumgsb_left%*%(t(cumgsb_all)*((delI*(eCHaz*eCHazLeft))/(exp2denfix^2)))
  I_part_6_ht = -cumgsb_all%*%(t(cumgsb_all)*((delI*(eCHaz^2))/(exp2denfix^2)))
  
  I_part_ht = I_part_1_ht + I_part_2_ht + I_part_3_ht + I_part_4_ht + I_part_5_ht + I_part_6_ht
  
  #penalty_dd_ht = 2*h*Rmat
  
  hess_theta = E_part_ht + R_part_ht + L_part_ht + I_part_ht  #- penalty_dd_ht
  
  ##############################
  ## Hessian of Beta & Gamma: ##
  ##############################
  
  kp_right = as.vector(kp_right)
  kp_d_right = as.vector(kp_d_right)
  kp_2d_right = as.vector(kp_2d_right)
  
  kp_left_int = as.vector(kp_left_int)
  kp_d_left_int = as.vector(kp_d_left_int)
  kp_2d_left_int = as.vector(kp_2d_left_int)
  
  E_part_hbg = delE*((-haz*haz_d - haz*haz_2d*kp_right 
                      + (haz_d^2)*kp_right)/(haz^2) + haz + haz_d*kp_right)*kp_d_right
  
  R_part_hbg = delR*(haz_d*kp_right + haz)*kp_d_right
  
  L_part_hbg = -delL*((eCHaz*(haz + haz_d*kp_right))/exp1den - ((haz^2)*kp_right*eCHaz)/(exp1den^2))*kp_d_right
  
  I_part_1_hbg = ((eCHazLeft*(kp_left_int*(haz_d_left - haz_left^2) + haz_left))/exp2denfix)*kp_d_left_int
  I_part_2_hbg = ((eCHaz*(kp_right*(haz_d - haz^2) + haz))/exp2denfix)*kp_d_right
  
  com_hbg = eCHaz*haz*kp_right - eCHazLeft*haz_left*kp_left_int
  
  I_part_3_hbg = ((com_hbg*haz*eCHaz)/(exp2denfix^2))*kp_d_right
  I_part_4_hbg = ((com_hbg*haz_left*eCHazLeft)/(exp2denfix^2))*kp_d_left_int
  
  I_part_hbg = delI*(I_part_1_hbg - I_part_2_hbg + I_part_3_hbg - I_part_4_hbg)
  
  hess_beta_gamma = t(XmatHess)%*%(E_part_hbg + R_part_hbg + L_part_hbg + I_part_hbg)
  
  ##############################
  ## Hessian of Beta & Theta: ##
  ##############################
  
  E_part_hbt = -delE*((haz*t(gsb_d_all) - haz_d*t(gsb_all))/(haz^2) - t(gsb_all))*kp_right
  R_part_hbt = delR*kp_right*t(gsb_left)
  L_part_hbt = -delL*((eCHaz*(t(gsb_all) - t(cumgsb_all)*haz))/exp1den - (eCHaz^2*haz*t(cumgsb_all))/(exp1den^2))*kp_right
  I_part_hbt = delI*((eCHazLeft*(t(gsb_left) - t(cumgsb_left)*haz_left)*kp_left_int)/exp2denfix
                 - (eCHaz*(t(gsb_all) - t(cumgsb_all)*haz)*kp_right)/exp2denfix
                 + ((eCHaz*haz*kp_right - eCHazLeft*haz_left*kp_left_int)*(t(cumgsb_all)*eCHaz- t(cumgsb_left)*eCHazLeft))/(exp2denfix^2)
  )
  
  hess_beta_theta = t(XmatHess)%*%(E_part_hbt + R_part_hbt + L_part_hbt + I_part_hbt)
  
  ###############################
  ## Hessian of Gamma & Theta: ##
  ###############################
  
  E_part_hgt = t(kp_d_right)%*%(delE*((haz*t(gsb_d_all) - haz_d*t(gsb_all))/(haz^2) - t(gsb_all)))
  R_part_hgt = -t(kp_d_right)%*%(delR*t(gsb_left))
  L_part_hgt = t(kp_d_right)%*%(delL*((eCHaz*(t(gsb_all) - t(cumgsb_all)*haz))/exp1den - (eCHaz^2*haz*t(cumgsb_all))/(exp1den^2)))
  
  I_part_1_hgt = t(kp_d_right)%*%(delI*(eCHaz*(t(gsb_all) - t(cumgsb_all)*haz))/exp2denfix)
  I_part_2_hgt = t(kp_d_left_int)%*%(-delI*(eCHazLeft*(t(gsb_left) - t(cumgsb_left)*haz_left))/exp2denfix)
  
  com2_hgt = t(cumgsb_all)*eCHaz - t(cumgsb_left)*eCHazLeft
  
  I_part_3_hgt = t(kp_d_left_int)%*%((delI*eCHazLeft*haz_left*com2_hgt)/(exp2denfix^2))
  I_part_4_hgt = t(kp_d_right)%*%((-delI*eCHaz*haz*com2_hgt)/(exp2denfix^2))
  
  I_part_hgt = I_part_1_hgt + I_part_2_hgt + I_part_3_hgt + I_part_4_hgt
  
  hess_gamma_theta = E_part_hgt + R_part_hgt + L_part_hgt + I_part_hgt
  
  #########################
  ## Final Full Hessian: ##
  #########################
  
  row1 = cbind(hess_beta, hess_beta_gamma, hess_beta_theta)
  row2 = cbind(t(hess_beta_gamma), hess_gamma, hess_gamma_theta)
  row3 = cbind(t(hess_beta_theta), t(hess_gamma_theta), hess_theta)
  
  fullHess = as.matrix(rbind(row1, row2, row3))
  rownames(fullHess) <- NULL
  colnames(fullHess) <- NULL
  
  return(fullHess = fullHess)
  
}