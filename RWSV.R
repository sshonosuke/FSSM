library(MCMCpack)
library(MASS)
library(isotone)
library(tidyverse)

run_rwsv_dlm <- function(Y, delta = 0.85, beta = 0.95, nsim = 3000, seed = 1234) {
  
  set.seed(seed)
  
  # --- Model Setup ---
  T <- nrow(Y)
  K <- ncol(Y)
  x_points <- seq(1/(1+K), 1-1/(1+K), length.out = K)
  
  if(K > 1/(1-beta)){ beta <- 0.98 }
  
  # Priors 
  m0 <- rep(0, K)
  c0 <- 10
  n0 <- (3 + K - 1)/2
  h0 <- n0 + K - 1
  D0 <- 0.01 * diag(K)
  
  # --- Internal Helpers ---
  
  # Wishart sampler for low degrees of freedom
  mywishart <- function(nu, Sigma){
    dim_Sigma <- nrow(Sigma)
    if (nu < dim_Sigma) {
      # Fallback for singular cases
      tryCatch({
        x <- mvrnorm(n = max(1, round(nu)), mu = rep(0, dim_Sigma), Sigma = Sigma)
        if(is.vector(x)) return(x %*% t(x)) else return(t(x) %*% x)
      }, error = function(e){
        warning("Wishart sampling failed, returning scaled identity.")
        return(diag(dim_Sigma) * 1e-6)
      })
    } else {
      return(rwish(nu, Sigma))
    }
  }
  
  # Trapezoidal Gini calculation
  calc_gini <- function(vals, x_pts){
    y_aug <- c(0, vals, 1)
    x_aug <- c(0, x_pts, 1)
    area <- sum(diff(x_aug) * (head(y_aug, -1) + tail(y_aug, -1)) / 2)
    return(1 - 2 * area)
  }
  
  # --- Forward Filtering ---
  # Storage
  m_list <- array(0, dim=c(K, 1, T))
  c_list <- rep(0, T)
  D_list <- array(0, dim=c(K, K, T))
  h_list <- rep(0, T) # degrees of freedom
  e_list <- array(0, dim=c(K, 1, T)) # forecast errors
  r_list <- rep(0, T)
  q_list <- rep(0, T)
  A_list <- rep(0, T) # adaptive coefficient
  
  # Initialization at t=1
  e_list[,,1] <- matrix(Y[1,] - m0, ncol=1)
  r_list[1]   <- c0 / delta
  q_list[1]   <- r_list[1] + 1
  A_list[1]   <- r_list[1] / q_list[1]
  D_list[,,1] <- beta * D0 + (e_list[,,1] %*% t(e_list[,,1])) / q_list[1]
  h_list[1]   <- beta * h0 + 1
  m_list[,,1] <- m0 + A_list[1] * e_list[,,1] 
  c_list[1]   <- r_list[1] - A_list[1]^2 * q_list[1]
  
  # Loop t=2 to T
  for(t in 2:T){
    e_list[,,t] <- matrix(Y[t,] - m_list[,,t-1], ncol=1)
    r_list[t]   <- c_list[t-1] / delta
    q_list[t]   <- r_list[t] + 1
    A_list[t]   <- r_list[t] / q_list[t]
    D_list[,,t] <- beta * D_list[,,t-1] + (e_list[,,t] %*% t(e_list[,,t])) / q_list[t]
    h_list[t]   <- beta * h_list[t-1] + 1
    m_list[,,t] <- m_list[,,t-1] + A_list[t] * e_list[,,t] 
    c_list[t]   <- r_list[t] - A_list[t]^2 * q_list[t]
  }
  
  # --- Backward Smoothing ---
  aT_list <- array(0, dim=c(K, 1, T))
  rT_list <- rep(0, T)
  
  aT_list[,,T] <- m_list[,,T]
  rT_list[T]   <- c_list[T]
  
  # Iterative smoothing
  check <- 1
  iter_check <- 0
  max_check <- 10
  
  while(check >= 1 && iter_check < max_check){
    for(t in (T-1):1){
      aT_list[,,t] <- (1-delta)*m_list[,,t] + delta*aT_list[,,t+1]
      rT_list[t]   <- (1-delta)*c_list[t]   + delta^2 * rT_list[t+1]
    }
    check <- sum(aT_list < 0) # Check
    iter_check <- iter_check + 1
  }
  
  # --- PAVA ---
  dlm_mean_projected <- matrix(NA, nrow=T, ncol=K)
  
  for(t in 1:T){
    pava_res <- gpava(c(0, x_points, 1), c(0, aT_list[,,t], 1))
    
    # Extract corrected inner points
    corrected <- matrix(pava_res$x[2:(K+1)], ncol=1)
    aT_list[,,t] <- corrected
    dlm_mean_projected[t,] <- corrected
  }
  
  # --- Retrospective Simulation ---
  Theta_list <- array(0, dim=c(K, 1, T, nsim))
  Sigma_list <- array(0, dim=c(K, K, T, nsim))
  y_rep      <- array(NA, dim=c(nsim, T, K)) # For predictive checks
  
  mstar_list <- array(0, dim=c(K, 1, T))
  cstar_list <- rep(0, T)
  
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  
  for(i in 1:nsim){
    
    # Initialization at T
    Sigma_list[,,T,i] <- riwish(h_list[T], D_list[,,T])
    Theta_list[,,T,i] <- mvrnorm(1, m_list[,,T], c_list[T] * Sigma_list[,,T,i])
    y_rep[i,T,]       <- mvrnorm(1, Theta_list[,,T,i], Sigma_list[,,T,i])
    
    mstar_list[,,T] <- m_list[,,T]
    cstar_list[T]   <- c_list[T]
    
    # Backward recursion
    for(t in (T-1):1){
      # Update moments
      mstar_list[,,t] <- m_list[,,t] + delta * (Theta_list[,,t+1,i] - m_list[,,t])
      cstar_list[t]   <- c_list[t] - delta^2 * r_list[t+1] 
      
      # Sample Volatility
      nu_val <- max(2, round((1-beta)*h_list[t])) 
      S_val  <- solve(D_list[,,t])
      S_val <- (S_val + t(S_val)) / 2 
      Gam    <- mywishart(nu_val, S_val)
      Sigma_list[,,t,i] <- solve( beta * solve(Sigma_list[,,t+1,i]) + Gam )
      
      # Sample State
      Theta_list[,,t,i] <- mvrnorm(1, mstar_list[,,t], cstar_list[t] * Sigma_list[,,t+1,i])
      
      # Sample Predictive
      y_rep[i,t,] <- mvrnorm(1, Theta_list[,,t,i], Sigma_list[,,t,i])
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # --- Post-Processing & Metrics ---
  
  # Gini Calculation (Trapezoidal on Projected Mean)
  gini_dlm <- numeric(T)
  for(t in 1:T){
    gini_dlm[t] <- calc_gini(dlm_mean_projected[t,], x_points)
  }
  
  # Log Predictive Power Score (LogPPL)
  logPPL_0 <- 0; logPPL_1 <- 0; logPPL_inf <- 0
  y_rep_mean <- apply(y_rep, c(2,3), mean)
  
  for(t in 1:T){
    for(k in 1:K){
      pred_var <- mean((y_rep[,t,k] - y_rep_mean[t,k])^2)
      bias_sq  <- (Y[t,k] - y_rep_mean[t,k])^2
      
      logPPL_0   <- logPPL_0 + pred_var
      logPPL_1   <- logPPL_1 + pred_var + 0.5 * bias_sq
      logPPL_inf <- logPPL_inf + pred_var + bias_sq
    }
  }
  logPPL_vals <- c(r0 = log(logPPL_0), r1 = log(logPPL_1), rInf = log(logPPL_inf))
  
  return(list(
    mean_projected = dlm_mean_projected,
    theta_samples  = Theta_list, # Note: These are unconstrained samples
    sigma_samples  = Sigma_list,
    gini           = gini_dlm,
    log_ppl        = logPPL_vals,
    Y_actual       = Y
  ))
}





# Example Usage
if (FALSE) {
  # Generate dummy data
  T_dummy <- 100
  raw_mat <- matrix(runif(T_dummy * 5), nrow = T_dummy, ncol = 5)
  raw_mat <- t(apply(raw_mat, 1, function(x) x / sum(x))) 
  cum_data <- t(apply(raw_mat, 1, cumsum))
  Y_dummy <- cum_data[, 1:4] 
  
  # Run DLM
  res <- run_rwsv_dlm(Y_dummy, delta = 0.85, beta = 0.95, nsim = 500)
  
  # Plot Gini
  plot(res$gini, type='l', main="Estimated Gini Coefficient", ylab="Gini")
  
  # Compare Projected Mean vs Observed at t=50
  t_check <- 50
  x_pts <- seq(0.2, 0.8, length.out=4)
  plot(x_pts, Y_dummy[t_check,], ylim=c(0,1), pch=19, main="Fit Check")
  lines(x_pts, res$mean_projected[t_check,], col="red", lwd=2)
}