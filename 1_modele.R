###############################################################################
################### 1 : MODEL DEFINITION - NEW ###############################
###############################################################################

###############################################################################
# ---- HELPER FUNCTIONS FOR THE MODEL ----
###############################################################################

# Compute total population counts for a given setting (hospital or community)
compute_totals <- function(S0, SA, S_II, S_III, C0, CA, C_II, C_III, I, I_II, I_III) {
  totals <- list(S = S0 + SA + S_II + S_III,
                 C = C0 + CA + C_II + C_III,
                 I_tot = I + I_II + I_III,
                 N = S0 + SA + S_II + S_III + C0 + CA + C_II + C_III + I + I_II + I_III)
  return(totals)
}

# Compute the force of infection (lambda)
compute_lambda <- function(beta, C, I_tot, N, nu) {
  lambda <- beta * (C + nu * I_tot) / N
  return(lambda)
}

# Compute stage-specific progression rates from sigmas
compute_sigmas <- function(sigma_base, k_A, k_II, k_III) {
  sigmas <- list(sigma_A = k_A * sigma_base,
                 sigma_II = k_II * sigma_base,
                 sigma_III = k_III * sigma_base)
  return(sigmas)
}

# Compute dynamic alpha rates at a given model state
compute_alpha_dynamic <- function(tot_h, tot_c, I_c, I_II_c, I_III_c,
                                  delta, w, w_I, w_II, w_III) {
  
  out_hc <- delta * (tot_h$S + tot_h$C) # hospital outflow, excluding infected compartments
  den_alpha <- w * (tot_c$S + tot_c$C) + w_I * I_c + w_II * I_II_c + w_III * I_III_c # normalization (with weights)

  if (den_alpha <= 0) den_alpha <- 1 # safety
  alpha <- out_hc / den_alpha
  
  alpha_rates <- list(alpha = alpha,
                      alpha_I = alpha * w_I,
                      alpha_II = alpha * w_II,
                      alpha_III = alpha * w_III)
  return(alpha_rates)
}

###############################################################################
# ---- CLOSTRIDIUM DIFFICILE TRANSMISSION MODEL ----
###############################################################################

cdiff_model <- function(t, pop, params) {
  
  # Default behavior: keep dynamic alpha unless explicitly switched to fixed
  use_fixed_alpha <- FALSE
  # if the parameter vector is named and contains "use_fixed_alpha" -> read it and activate fixed alpha only when its value is exactly 1
  if (!is.null(names(params)) && ("use_fixed_alpha" %in% names(params))) {
    use_fixed_alpha <- (as.numeric(params[["use_fixed_alpha"]]) == 1)
  }

  derivatives <- with(as.list(c(pop, params)), {
    
    # Totals hospital / community 
    tot_h <- compute_totals(S0_h, SA_h, S_II_h, S_III_h, C0_h, CA_h, C_II_h, C_III_h, I_h, I_II_h, I_III_h)
    tot_c <- compute_totals(S0_c, SA_c, S_II_c, S_III_c, C0_c, CA_c, C_II_c, C_III_c, I_c, I_II_c, I_III_c)
    
    # Alpha choice:
    # - default: dynamic alpha computed at each time step -> FOR CALIBRATION 
    # - optional: fixed alpha from params if use_fixed_alpha = 1 -> FOR COMPARAISON SIMULATION
    if (use_fixed_alpha) {
      alpha <- alpha_const
      alpha_I <- alpha * w_I
      alpha_II <- alpha * w_II
      alpha_III <- alpha * w_III
    } else {
      alpha_rates <- compute_alpha_dynamic(tot_h = tot_h, tot_c = tot_c, I_c = I_c, I_II_c = I_II_c, I_III_c = I_III_c,
        delta = delta, w = w, w_I = w_I, w_II = w_II, w_III = w_III)
      alpha <- alpha_rates$alpha
      alpha_I <- alpha_rates$alpha_I
      alpha_II <- alpha_rates$alpha_II
      alpha_III <- alpha_rates$alpha_III
    }
    
    # Forces of infection 
    lambda_h <- compute_lambda(beta_h, tot_h$C, tot_h$I_tot, tot_h$N, nu)
    lambda_c <- compute_lambda(beta_c, tot_c$C, tot_c$I_tot, tot_c$N, nu)
    
    # Sigmas 
    sig_h <- compute_sigmas(sigma_h, k_A, k_II, k_III)
    sig_c <- compute_sigmas(sigma_c, k_A, k_II, k_III)
    
    # ---- EDO - Hospital ----
    
    # Primary infection
    dS0_h <- -lambda_h*S0_h + gamma*C0_h - tau_h*S0_h + omega*SA_h + phi*(S_II_h + S_III_h) + alpha*S0_c - delta*S0_h
    
    dSA_h <- -lambda_h*SA_h + gamma*CA_h + tau_h*S0_h - omega*SA_h +  alpha*SA_c - delta*SA_h
    
    dC0_h <- lambda_h*S0_h - gamma*C0_h - tau_h*C0_h + omega*CA_h - sigma_h*C0_h + alpha*C0_c - delta*C0_h
    
    dCA_h <- lambda_h*SA_h - gamma*CA_h + tau_h*C0_h - omega*CA_h - sig_h$sigma_A*CA_h + alpha*CA_c - delta*CA_h
    
    dI_h <- sigma_h*C0_h + sig_h$sigma_A*CA_h - epsilon*I_h + alpha_I*I_c
    
    # Second episode
    dS_II_h <- p*epsilon*I_h - lambda_h*S_II_h + gamma*C_II_h - phi*S_II_h + alpha*S_II_c - delta*S_II_h
    
    dC_II_h <- (1-p)*epsilon*I_h + lambda_h*S_II_h - gamma*C_II_h - sig_h$sigma_II*C_II_h + alpha*C_II_c - delta*C_II_h
    
    dI_II_h <- sig_h$sigma_II*C_II_h - epsilon*I_II_h + alpha_II*I_II_c
    
    # Third or more episode
    dS_III_h <- p*epsilon*(I_II_h + I_III_h) - lambda_h*S_III_h + gamma*C_III_h - phi*S_III_h + alpha*S_III_c - delta*S_III_h
    
    dC_III_h <- (1-p)*epsilon*(I_II_h + I_III_h) + lambda_h*S_III_h - gamma*C_III_h - sig_h$sigma_III*C_III_h + alpha*C_III_c - delta*C_III_h
    
    dI_III_h <- sig_h$sigma_III*C_III_h - epsilon*I_III_h + alpha_III*I_III_c
    
    # ---- EDO - Community ----
    
    # Primary infection
    dS0_c <- -lambda_c*S0_c + gamma*C0_c - tau_c*S0_c + omega*SA_c + phi*(S_II_c + S_III_c) - alpha*S0_c + delta*S0_h
    
    dSA_c <- -lambda_c*SA_c + gamma*CA_c + tau_c*S0_c - omega*SA_c - alpha*SA_c + delta*SA_h
    
    dC0_c <- lambda_c*S0_c - gamma*C0_c - tau_c*C0_c + omega*CA_c - sigma_c*C0_c - alpha*C0_c + delta*C0_h
    
    dCA_c <- lambda_c*SA_c - gamma*CA_c + tau_c*C0_c - omega*CA_c - sig_c$sigma_A*CA_c - alpha*CA_c + delta*CA_h
    
    dI_c <- sigma_c*C0_c + sig_c$sigma_A*CA_c - epsilon*I_c - alpha_I*I_c 
    
    # Second episode (first recurrence)
    dS_II_c <- p*epsilon*I_c - lambda_c*S_II_c + gamma*C_II_c - phi*S_II_c - alpha*S_II_c + delta*S_II_h
    
    dC_II_c <- (1-p)*epsilon*I_c + lambda_c*S_II_c - gamma*C_II_c - sig_c$sigma_II*C_II_c - alpha*C_II_c + delta*C_II_h
    
    dI_II_c <- sig_c$sigma_II*C_II_c - epsilon*I_II_c - alpha_II*I_II_c 
    
    # Third+ episode (second+ recurrence)
    dS_III_c <- p*epsilon*(I_II_c + I_III_c) - lambda_c*S_III_c  + gamma*C_III_c - phi*S_III_c - alpha*S_III_c + delta*S_III_h
    
    dC_III_c <- (1-p)*epsilon*(I_II_c + I_III_c) + lambda_c*S_III_c - gamma*C_III_c - sig_c$sigma_III*C_III_c - alpha*C_III_c + delta*C_III_h
    
    dI_III_c <- sig_c$sigma_III*C_III_c - epsilon*I_III_c - alpha_III*I_III_c 

    list(c(
      S0_h=dS0_h, SA_h=dSA_h, C0_h=dC0_h, CA_h=dCA_h, I_h=dI_h,
      S_II_h=dS_II_h, C_II_h=dC_II_h, I_II_h=dI_II_h,
      S_III_h=dS_III_h, C_III_h=dC_III_h, I_III_h=dI_III_h,
      S0_c=dS0_c, SA_c=dSA_c, C0_c=dC0_c, CA_c=dCA_c, I_c=dI_c,
      S_II_c=dS_II_c, C_II_c=dC_II_c, I_II_c=dI_II_c,
      S_III_c=dS_III_c, C_III_c=dC_III_c, I_III_c=dI_III_c))
  })
  return(derivatives)
}


###############################################################################
# ---- VACCINATION MODEL (stratified by vaccination status) ----
###############################################################################

cdiff_vacc <- function(t, pop, params) {
  derivatives <- with(as.list(c(pop, params)), {
    
    # Totals for non-vaccinated (nv)
    tot_h_nv <- compute_totals(S0_h_nv, SA_h_nv, S_II_h_nv, S_III_h_nv, C0_h_nv, CA_h_nv, C_II_h_nv, C_III_h_nv, I_h_nv, I_II_h_nv, I_III_h_nv)
    tot_c_nv <- compute_totals(S0_c_nv, SA_c_nv, S_II_c_nv, S_III_c_nv, C0_c_nv, CA_c_nv, C_II_c_nv, C_III_c_nv, I_c_nv, I_II_c_nv, I_III_c_nv)
    # Totals for vaccinated (v)
    tot_h_v <- compute_totals(S0_h_v, SA_h_v, S_II_h_v, S_III_h_v, C0_h_v, CA_h_v, C_II_h_v, C_III_h_v, I_h_v, I_II_h_v, I_III_h_v)
    tot_c_v <- compute_totals(S0_c_v, SA_c_v, S_II_c_v, S_III_c_v, C0_c_v, CA_c_v, C_II_c_v, C_III_c_v, I_c_v, I_II_c_v, I_III_c_v)
    
    # Fixed alpha values (same for all vaccination scenarios)
    alpha <- alpha_const
    alpha_I <- alpha * w_I
    alpha_II <- alpha * w_II
    alpha_III <- alpha * w_III
    
    # Common force of infection (FOI) 
    C_h_total <- tot_h_nv$C + tot_h_v$C
    I_h_total <- tot_h_nv$I_tot + tot_h_v$I_tot
    N_h_total <- tot_h_nv$N + tot_h_v$N
    lambda_h <- compute_lambda(beta_h, C_h_total, I_h_total, N_h_total, nu)
    
    C_c_total <- tot_c_nv$C + tot_c_v$C
    I_c_total <- tot_c_nv$I_tot + tot_c_v$I_tot
    N_c_total <- tot_c_nv$N + tot_c_v$N
    lambda_c <- compute_lambda(beta_c, C_c_total, I_c_total, N_c_total, nu)
    
    # Sigmas (with vaccine effect for vaccinated)
    sig_h_nv <- compute_sigmas(sigma_h, k_A, k_II, k_III)
    sig_c_nv <- compute_sigmas(sigma_c, k_A, k_II, k_III)
    
    sigma_h_v <- sigma_h * (1 - VE)
    sigma_c_v <- sigma_c * (1 - VE)
    sig_h_v <- compute_sigmas(sigma_h_v, k_A, k_II, k_III)
    sig_c_v <- compute_sigmas(sigma_c_v, k_A, k_II, k_III)
    
    # ---- EDO - NON-VACCINATED POPULATION ----
    
    # Hospital - non-vaccinated
    dS0_h_nv <- -lambda_h*S0_h_nv + gamma*C0_h_nv - tau_h*S0_h_nv + omega*SA_h_nv + phi*(S_II_h_nv + S_III_h_nv) + alpha*S0_c_nv - delta*S0_h_nv
    
    dSA_h_nv <- -lambda_h*SA_h_nv + gamma*CA_h_nv + tau_h*S0_h_nv - omega*SA_h_nv + alpha*SA_c_nv - delta*SA_h_nv
    
    dC0_h_nv <- lambda_h*S0_h_nv - gamma*C0_h_nv - tau_h*C0_h_nv + omega*CA_h_nv - sigma_h*C0_h_nv + alpha*C0_c_nv - delta*C0_h_nv
    
    dCA_h_nv <- lambda_h*SA_h_nv - gamma*CA_h_nv + tau_h*C0_h_nv - omega*CA_h_nv - sig_h_nv$sigma_A*CA_h_nv + alpha*CA_c_nv - delta*CA_h_nv
    
    dI_h_nv <- sigma_h*C0_h_nv + sig_h_nv$sigma_A*CA_h_nv - epsilon*I_h_nv + alpha_I*I_c_nv
    
    dS_II_h_nv <- p*epsilon*I_h_nv + gamma*C_II_h_nv - lambda_h*S_II_h_nv - phi*S_II_h_nv + alpha*S_II_c_nv - delta*S_II_h_nv
    
    dC_II_h_nv <- (1-p)*epsilon*I_h_nv + lambda_h*S_II_h_nv - gamma*C_II_h_nv - sig_h_nv$sigma_II*C_II_h_nv + alpha*C_II_c_nv - delta*C_II_h_nv
    
    dI_II_h_nv <- sig_h_nv$sigma_II*C_II_h_nv - epsilon*I_II_h_nv + alpha_II*I_II_c_nv
    
    dS_III_h_nv <- p*epsilon*(I_II_h_nv + I_III_h_nv) + gamma*C_III_h_nv - lambda_h*S_III_h_nv - phi*S_III_h_nv + alpha*S_III_c_nv - delta*S_III_h_nv
    
    dC_III_h_nv <- (1-p)*epsilon*(I_II_h_nv + I_III_h_nv) + lambda_h*S_III_h_nv - gamma*C_III_h_nv - sig_h_nv$sigma_III*C_III_h_nv + alpha*C_III_c_nv - delta*C_III_h_nv
    
    dI_III_h_nv <- sig_h_nv$sigma_III*C_III_h_nv - epsilon*I_III_h_nv + alpha_III*I_III_c_nv
    
    # Community - non-vaccinated
    dS0_c_nv <- -lambda_c*S0_c_nv + gamma*C0_c_nv - tau_c*S0_c_nv + omega*SA_c_nv + phi*(S_II_c_nv + S_III_c_nv) - alpha*S0_c_nv + delta*S0_h_nv
    
    dSA_c_nv <- -lambda_c*SA_c_nv + gamma*CA_c_nv + tau_c*S0_c_nv - omega*SA_c_nv - alpha*SA_c_nv + delta*SA_h_nv
    
    dC0_c_nv <- lambda_c*S0_c_nv - gamma*C0_c_nv - tau_c*C0_c_nv + omega*CA_c_nv - sigma_c*C0_c_nv - alpha*C0_c_nv + delta*C0_h_nv
    
    dCA_c_nv <- lambda_c*SA_c_nv - gamma*CA_c_nv + tau_c*C0_c_nv - omega*CA_c_nv - sig_c_nv$sigma_A*CA_c_nv - alpha*CA_c_nv + delta*CA_h_nv
    
    dI_c_nv <- sigma_c*C0_c_nv + sig_c_nv$sigma_A*CA_c_nv - epsilon*I_c_nv - alpha_I*I_c_nv
    
    dS_II_c_nv <- p*epsilon*I_c_nv + gamma*C_II_c_nv - lambda_c*S_II_c_nv - phi*S_II_c_nv - alpha*S_II_c_nv + delta*S_II_h_nv
    
    dC_II_c_nv <- (1-p)*epsilon*I_c_nv + lambda_c*S_II_c_nv - gamma*C_II_c_nv - sig_c_nv$sigma_II*C_II_c_nv - alpha*C_II_c_nv + delta*C_II_h_nv
    
    dI_II_c_nv <- sig_c_nv$sigma_II*C_II_c_nv - epsilon*I_II_c_nv - alpha_II*I_II_c_nv
    
    dS_III_c_nv <- p*epsilon*(I_II_c_nv + I_III_c_nv) + gamma*C_III_c_nv - lambda_c*S_III_c_nv - phi*S_III_c_nv - alpha*S_III_c_nv + delta*S_III_h_nv
    
    dC_III_c_nv <- (1-p)*epsilon*(I_II_c_nv + I_III_c_nv) + lambda_c*S_III_c_nv - gamma*C_III_c_nv - sig_c_nv$sigma_III*C_III_c_nv - alpha*C_III_c_nv + delta*C_III_h_nv
    
    dI_III_c_nv <- sig_c_nv$sigma_III*C_III_c_nv - epsilon*I_III_c_nv - alpha_III*I_III_c_nv
    
    # ---- EDO - VACCINATED POPULATION ----
    
    # Hospital - vaccinated
    dS0_h_v <- -lambda_h*S0_h_v + gamma*C0_h_v - tau_h*S0_h_v + omega*SA_h_v + phi*(S_II_h_v + S_III_h_v) + alpha*S0_c_v - delta*S0_h_v
    
    dSA_h_v <- -lambda_h*SA_h_v + gamma*CA_h_v + tau_h*S0_h_v - omega*SA_h_v + alpha*SA_c_v - delta*SA_h_v
    
    dC0_h_v <- lambda_h*S0_h_v - gamma*C0_h_v - tau_h*C0_h_v + omega*CA_h_v - sigma_h_v*C0_h_v + alpha*C0_c_v - delta*C0_h_v
    
    dCA_h_v <- lambda_h*SA_h_v - gamma*CA_h_v + tau_h*C0_h_v - omega*CA_h_v - sig_h_v$sigma_A*CA_h_v + alpha*CA_c_v - delta*CA_h_v
    
    dI_h_v <- sigma_h_v*C0_h_v + sig_h_v$sigma_A*CA_h_v - epsilon*I_h_v + alpha_I*I_c_v
    
    dS_II_h_v <- p*epsilon*I_h_v + gamma*C_II_h_v - lambda_h*S_II_h_v - phi*S_II_h_v + alpha*S_II_c_v - delta*S_II_h_v
    
    dC_II_h_v <- (1-p)*epsilon*I_h_v + lambda_h*S_II_h_v - gamma*C_II_h_v - sig_h_v$sigma_II*C_II_h_v + alpha*C_II_c_v - delta*C_II_h_v
    
    dI_II_h_v <- sig_h_v$sigma_II*C_II_h_v - epsilon*I_II_h_v + alpha_II*I_II_c_v
    
    dS_III_h_v <- p*epsilon*(I_II_h_v + I_III_h_v) + gamma*C_III_h_v - lambda_h*S_III_h_v - phi*S_III_h_v + alpha*S_III_c_v - delta*S_III_h_v
    
    dC_III_h_v <- (1-p)*epsilon*(I_II_h_v + I_III_h_v) + lambda_h*S_III_h_v - gamma*C_III_h_v - sig_h_v$sigma_III*C_III_h_v + alpha*C_III_c_v - delta*C_III_h_v
    
    dI_III_h_v <- sig_h_v$sigma_III*C_III_h_v - epsilon*I_III_h_v + alpha_III*I_III_c_v
    
    # Community - vaccinated
    dS0_c_v <- -lambda_c*S0_c_v + gamma*C0_c_v - tau_c*S0_c_v + omega*SA_c_v + phi*(S_II_c_v + S_III_c_v) - alpha*S0_c_v + delta*S0_h_v
    
    dSA_c_v <- -lambda_c*SA_c_v + gamma*CA_c_v + tau_c*S0_c_v - omega*SA_c_v - alpha*SA_c_v + delta*SA_h_v
    
    dC0_c_v <- lambda_c*S0_c_v - gamma*C0_c_v - tau_c*C0_c_v + omega*CA_c_v - sigma_c_v*C0_c_v - alpha*C0_c_v + delta*C0_h_v
    
    dCA_c_v <- lambda_c*SA_c_v - gamma*CA_c_v + tau_c*C0_c_v - omega*CA_c_v - sig_c_v$sigma_A*CA_c_v - alpha*CA_c_v + delta*CA_h_v
    
    dI_c_v <- sigma_c_v*C0_c_v + sig_c_v$sigma_A*CA_c_v - epsilon*I_c_v - alpha_I*I_c_v
    
    dS_II_c_v <- p*epsilon*I_c_v + gamma*C_II_c_v - lambda_c*S_II_c_v - phi*S_II_c_v - alpha*S_II_c_v + delta*S_II_h_v
    
    dC_II_c_v <- (1-p)*epsilon*I_c_v + lambda_c*S_II_c_v - gamma*C_II_c_v - sig_c_v$sigma_II*C_II_c_v - alpha*C_II_c_v + delta*C_II_h_v
    
    dI_II_c_v <- sig_c_v$sigma_II*C_II_c_v - epsilon*I_II_c_v - alpha_II*I_II_c_v
    
    dS_III_c_v <- p*epsilon*(I_II_c_v + I_III_c_v) + gamma*C_III_c_v - lambda_c*S_III_c_v - phi*S_III_c_v - alpha*S_III_c_v + delta*S_III_h_v
    
    dC_III_c_v <- (1-p)*epsilon*(I_II_c_v + I_III_c_v) + lambda_c*S_III_c_v - gamma*C_III_c_v - sig_c_v$sigma_III*C_III_c_v - alpha*C_III_c_v + delta*C_III_h_v
    
    dI_III_c_v <- sig_c_v$sigma_III*C_III_c_v - epsilon*I_III_c_v - alpha_III*I_III_c_v
    d_state <- c(
      # Non-vaccinated
      S0_h_nv=dS0_h_nv, SA_h_nv=dSA_h_nv, C0_h_nv=dC0_h_nv, CA_h_nv=dCA_h_nv, I_h_nv=dI_h_nv,
      S_II_h_nv=dS_II_h_nv, C_II_h_nv=dC_II_h_nv, I_II_h_nv=dI_II_h_nv,
      S_III_h_nv=dS_III_h_nv, C_III_h_nv=dC_III_h_nv, I_III_h_nv=dI_III_h_nv,
      S0_c_nv=dS0_c_nv, SA_c_nv=dSA_c_nv, C0_c_nv=dC0_c_nv, CA_c_nv=dCA_c_nv, I_c_nv=dI_c_nv,
      S_II_c_nv=dS_II_c_nv, C_II_c_nv=dC_II_c_nv, I_II_c_nv=dI_II_c_nv,
      S_III_c_nv=dS_III_c_nv, C_III_c_nv=dC_III_c_nv, I_III_c_nv=dI_III_c_nv,
      # Vaccinated
      S0_h_v=dS0_h_v, SA_h_v=dSA_h_v, C0_h_v=dC0_h_v, CA_h_v=dCA_h_v, I_h_v=dI_h_v,
      S_II_h_v=dS_II_h_v, C_II_h_v=dC_II_h_v, I_II_h_v=dI_II_h_v,
      S_III_h_v=dS_III_h_v, C_III_h_v=dC_III_h_v, I_III_h_v=dI_III_h_v,
      S0_c_v=dS0_c_v, SA_c_v=dSA_c_v, C0_c_v=dC0_c_v, CA_c_v=dCA_c_v, I_c_v=dI_c_v,
      S_II_c_v=dS_II_c_v, C_II_c_v=dC_II_c_v, I_II_c_v=dI_II_c_v,
      S_III_c_v=dS_III_c_v, C_III_c_v=dC_III_c_v, I_III_c_v=dI_III_c_v
    )

    list(d_state)
  })
  return(derivatives)
}


###############################################################################
# ---- RUN MODEL EXECUTION ----
###############################################################################

# Run the ODE model over a fixed time grid (used inside run_model_until_equilibrium and for scenario simulations)
run_model <- function(params_vec, init_cond, time_vec, model_function = cdiff_model) {
  out <- as.data.frame(
    lsoda(
      y = init_cond,
      times = time_vec,
      func = model_function,
      parms = params_vec,
      atol = 1e-10,
      rtol = 1e-10,
      maxsteps = 100000
    )
  )
  return(out)
}

# Run the ODE model until equilibrium or until the maximum time is reached
run_model_until_equilibrium <- function(params_vec, init_cond, time_max, by = 1,
                                        equilibrium_tol = 1e-6,
                                        min_time_before_check = 365,
                                        chunk_length = 365,
                                        model_function = cdiff_model) {
  out_all <- NULL # stores the full simulation output
  current_state <- init_cond # current initial state for each chunk
  reached_equilibrium <- FALSE # becomes TRUE if equilibrium is reached
  last_max_abs_dydt <- NA_real_ # stores the last max |dX/dt|

  for (start_time in seq(0, time_max, by = chunk_length)) { # run chunk by chunk
    if (start_time >= time_max) break # stop if maximum time is reached

    end_time <- min(start_time + chunk_length, time_max) # end of current chunk
    times_chunk <- unique(c(seq(start_time, end_time, by = by), end_time)) # times for this chunk

    out_chunk <- run_model( # simulate this chunk
      params_vec = params_vec, # model parameters
      init_cond = current_state, # start from current state
      time_vec = times_chunk, # time points for this chunk
      model_function = model_function # model to run
    )

    if (is.null(out_all)) { # first chunk
      out_all <- out_chunk # save all rows
    } else {
      out_all <- rbind(out_all, out_chunk[-1, , drop = FALSE]) # add rows, without duplicate time
    }

    current_state <- as.numeric(out_chunk[nrow(out_chunk), -1]) # last simulated state
    names(current_state) <- colnames(out_chunk)[-1] # keep compartment names

    if (end_time >= min_time_before_check) { # wait before checking equilibrium
      dy_eq <- model_function(0, current_state, params_vec)[[1]] # derivatives at current state
      last_max_abs_dydt <- max(abs(dy_eq)) # equilibrium criterion

      if (last_max_abs_dydt <= equilibrium_tol) { # equilibrium reached
        reached_equilibrium <- TRUE # save diagnostic
        break # stop simulation
      }
    }
  }

  attr(out_all, "equilibrium_reached") <- reached_equilibrium # TRUE/FALSE
  attr(out_all, "equilibrium_tol") <- equilibrium_tol # threshold used
  attr(out_all, "equilibrium_max_abs_dydt") <- last_max_abs_dydt # final max |dX/dt|

  return(out_all) # return simulation plus diagnostics
}


###############################################################################
# ---- OUTPUT METRICS ----
###############################################################################

# Compute total population size in hospital, community, or both
compute_population_totals <- function(last_state, setting = c("h", "c", "both")) {
  setting <- match.arg(setting)

  N_h <- last_state$S0_h + last_state$SA_h + last_state$S_II_h + last_state$S_III_h + last_state$C0_h + last_state$CA_h + last_state$C_II_h + last_state$C_III_h + last_state$I_h + last_state$I_II_h + last_state$I_III_h
  N_c <- last_state$S0_c + last_state$SA_c + last_state$S_II_c + last_state$S_III_c + last_state$C0_c + last_state$CA_c + last_state$C_II_c + last_state$C_III_c + last_state$I_c + last_state$I_II_c + last_state$I_III_c

  if (setting == "h") return(N_h)
  if (setting == "c") return(N_c)
  if (setting == "both") return(N_h + N_c)
}

# Compute asymptomatic carriage prevalence by setting
compute_carriage_prevalence <- function(last_state, setting = c("h", "c", "both")) {
  setting <- match.arg(setting)

  C_h <- last_state$C0_h + last_state$CA_h + last_state$C_II_h + last_state$C_III_h
  C_c <- last_state$C0_c + last_state$CA_c + last_state$C_II_c + last_state$C_III_c

  N_h <- compute_population_totals(last_state, "h")
  N_c <- compute_population_totals(last_state, "c")

  if (setting == "h") return(C_h / N_h)
  if (setting == "c") return(C_c / N_c)
  if (setting == "both") return((C_h + C_c) / (N_h + N_c))
}

# Compute instantaneous CDI incidence
compute_CDI_incidence <- function(last_state, params_vec, setting = c("h", "c", "both"), type = c("total", "primo", "recurrent")) {
  setting <- match.arg(setting)
  type <- match.arg(type)

  k_A <- as.numeric(params_vec["k_A"])
  k_II <- as.numeric(params_vec["k_II"])
  k_III <- as.numeric(params_vec["k_III"])
  sigma_h <- as.numeric(params_vec["sigma_h"])
  sigma_c <- as.numeric(params_vec["sigma_c"])

  calc_incidence <- function(sigma, C0, CA, C_II, C_III) {
    primo <- sigma * C0 + k_A * sigma * CA
    recurrent <- k_II * sigma * C_II + k_III * sigma * C_III
    total <- primo + recurrent
    return(c(total = total, primo = primo, recurrent = recurrent))
  }

  inc_h <- calc_incidence(sigma_h, last_state$C0_h, last_state$CA_h, last_state$C_II_h, last_state$C_III_h)
  inc_c <- calc_incidence(sigma_c, last_state$C0_c, last_state$CA_c, last_state$C_II_c, last_state$C_III_c)

  if (setting == "h") return(as.numeric(inc_h[type]))
  if (setting == "c") return(as.numeric(inc_c[type]))
  if (setting == "both") return(as.numeric((inc_h + inc_c)[type]))
}

# Compute recurrence prevalence ratios
compute_recurrence_prevalence <- function(last_state, setting = c("h", "c", "both"), type = c("rec_1", "rec_2")) {
  setting <- match.arg(setting)
  type <- match.arg(type)

  calc_recurrence <- function(I, I_II, I_III) {
    return(c(rec_1 = I_II / I, rec_2 = I_III / I_II))
  }

  rec_h <- calc_recurrence(last_state$I_h, last_state$I_II_h, last_state$I_III_h)
  rec_c <- calc_recurrence(last_state$I_c, last_state$I_II_c, last_state$I_III_c)
  rec_both <- calc_recurrence(last_state$I_h + last_state$I_c,
                              last_state$I_II_h + last_state$I_II_c,
                              last_state$I_III_h + last_state$I_III_c)

  if (setting == "h") return(rec_h[type])
  if (setting == "c") return(rec_c[type])
  if (setting == "both") return(rec_both[type])
}

# Compute the main model outputs from the final model state
compute_all_metrics <- function(last_state, params_vec, ode_result = NULL, targets = NULL) {
  if (!is.list(last_state)) {
    last_state <- as.list(last_state[1, ])
  }

  N_h <- compute_population_totals(last_state, "h")
  N_c <- compute_population_totals(last_state, "c")
  N_tot <- N_h + N_c

  metrics <- list(
    portage_h = compute_carriage_prevalence(last_state, "h"),
    portage_c = compute_carriage_prevalence(last_state, "c"),
    portage_tot = compute_carriage_prevalence(last_state, "both"),
    inc_h = unname(compute_CDI_incidence(last_state, params_vec, "h", "total")) / N_tot,
    inc_c = unname(compute_CDI_incidence(last_state, params_vec, "c", "total")) / N_tot,
    inc_tot = unname(compute_CDI_incidence(last_state, params_vec, "both", "total")) / N_tot,
    inc_primo_h = unname(compute_CDI_incidence(last_state, params_vec, "h", "primo")) / N_tot,
    inc_primo_c = unname(compute_CDI_incidence(last_state, params_vec, "c", "primo")) / N_tot,
    inc_primo_tot = unname(compute_CDI_incidence(last_state, params_vec, "both", "primo")) / N_tot,
    inc_rec_h = unname(compute_CDI_incidence(last_state, params_vec, "h", "recurrent")) / N_tot,
    inc_rec_c = unname(compute_CDI_incidence(last_state, params_vec, "c", "recurrent")) / N_tot,
    inc_rec_tot = unname(compute_CDI_incidence(last_state, params_vec, "both", "recurrent")) / N_tot,
    recid_1_tot = compute_recurrence_prevalence(last_state, "both", "rec_1"),
    recid_1_h = compute_recurrence_prevalence(last_state, "h", "rec_1"),
    recid_1_c = compute_recurrence_prevalence(last_state, "c", "rec_1"),
    recid_2_tot = compute_recurrence_prevalence(last_state, "both", "rec_2"),
    recid_2_h = compute_recurrence_prevalence(last_state, "h", "rec_2"),
    recid_2_c = compute_recurrence_prevalence(last_state, "c", "rec_2"),
    rec_burden_h = last_state$I_II_h + last_state$I_III_h,
    rec_burden_c = last_state$I_II_c + last_state$I_III_c,
    rec_burden_tot = last_state$I_II_h + last_state$I_III_h + last_state$I_II_c + last_state$I_III_c,
    N_h = N_h,
    N_c = N_c
  )

  if (!is.null(targets)) {
    metrics$errors <- list(
      err_portage_h = (metrics$portage_h - targets$portage_h) / targets$portage_h,
      err_portage_c = (metrics$portage_c - targets$portage_c) / targets$portage_c,
      err_inc_h = (metrics$inc_h - targets$incidence_h) / targets$incidence_h,
      err_inc_c = (metrics$inc_c - targets$incidence_c) / targets$incidence_c,
      err_recid_1_tot = (metrics$recid_1_tot - targets$recid_1) / targets$recid_1,
      err_recid_2_tot = (metrics$recid_2_tot - targets$recid_2) / targets$recid_2
    )
  }

  return(metrics)
}
