library(BayesLogit)
library(tidyverse)
library(truncnorm)
library(mvtnorm)
library(MASS)
library(MCMCpack)

# Shape-Constrained Functional Time Series MCMC (Improved FFBS)
run_shape_constrained_fssm <- function(Y, Psi, burnin = 20000, nsim = 80000, seed = 4321) {
  
  set.seed(seed)
  T <- nrow(Y)
  K <- ncol(Y)
  L <- ncol(Psi)
  
  # --- Model Setup ---
  
  # Hyperparameters
  a0 <- 0.8; v0 <- 0.2          
  B_a0 <- 0.0; B_v0 <- 5.0      
  und_d <- 0.001; und_s <- 0.001 
  n0 <- 0.001; s0 <- 0.001      
  
  # Sample AR(1) coef A and Sigma
  sample_A_Sigma_diag <- function(A, B, Sigma, u0, U){
    Alpha  <- t(t(U) - B)
    alpha0 <- u0 - B
    new_A <- numeric(L-1); new_Sigma <- numeric(L-1)
    
    x0_mat <- rbind(alpha0, Alpha[1:(T-1), , drop=FALSE])
    x1_mat <- Alpha
    
    for (l in 1:(L-1)){
      x0 <- x0_mat[,l]; x1 <- x1_mat[,l]
      sig_sq <- Sigma[l,l]
      
      denom <- sum(x0^2) / sig_sq + 1/v0^2
      var_phi <- 1 / denom
      mean_phi <- var_phi * (sum(x0 * x1) / sig_sq + a0/v0^2)
      
      old_phi <- A[l,l]
      cand_phi <- rtruncnorm(1, a=-1, b=1, mean=mean_phi, sd=sqrt(var_phi))
      
      if(runif(1) < sqrt((1 - cand_phi^2) / (1 - old_phi^2))){
        new_A[l] <- cand_phi
      } else {
        new_A[l] <- old_phi
      }
      
      phi_curr <- new_A[l]
      shape_post <- (T + und_d + 1) / 2
      rate_post  <- (sum((x1 - phi_curr * x0)^2) + und_s + (1 - phi_curr^2) * alpha0[l]^2) / 2
      new_Sigma[l] <- 1 / rgamma(1, shape = shape_post, rate = rate_post)
    }
    return(list(A = diag(new_A), Sigma = diag(new_Sigma)))
  }
  
  # Sample B
  sample_B <- function(A, Sigma, u0, U){
    sampleB <- numeric(L-1)
    U_lag <- rbind(u0, U[1:(T-1), , drop=FALSE])
    for (l in 1:(L-1)){
      phi <- A[l,l]; sig <- Sigma[l,l]
      y_vec <- U[,l] - phi * U_lag[,l]
      x_val <- (1 - phi)
      
      var_post  <- 1 / (1 / B_v0^2 + T * x_val^2 / sig)
      mean_post <- var_post * ((sum(y_vec) * x_val / sig) + (B_a0 / B_v0^2))
      
      # Correction
      prec_u0 <- (1 - phi^2) / sig
      var_exact <- 1 / (1/var_post + prec_u0)
      mean_exact <- var_exact * (mean_post/var_post + u0[l] * prec_u0)
      sampleB[l] <- rnorm(1, mean_exact, sqrt(var_exact))
    }
    return(sampleB)
  }
  
  # Sample sigma
  sample_sigma <- function(U){
    exp_U <- exp(U); denom <- 1 + rowSums(exp_U)
    Pi <- cbind(1/denom, exp_U/denom)
    ssr <- sum((Y - Pi %*% t(Psi))^2)
    return(sqrt(1 / rgamma(1, (n0 + T * K) / 2, (s0 + ssr) / 2)))
  }
  
  # Sample U & u0
  sample_U_scalar_FFBS <- function(U, u0, A, B, Sigma, sigma){
    new_u0 <- numeric(L-1)
    for (l in 1:(L-1)) {
      phi <- A[l, l]; sig_state <- Sigma[l, l]; mu <- B[l]
      ar_drift <- mu * (1 - phi)
      
      pseudo_y <- numeric(T); omega <- numeric(T)
      
      U_safe <- pmin(pmax(U, -25), 25)
      
      # -- PG Augmentation --
      for (t in 1:T) {
        idx <- l + 1
        D <- matrix(Y[t,], nrow=K, ncol=L) - Psi
        
        S_val <- 1 + sum(exp(U_safe[t, -l]))
        v_val <- exp(U_safe[t, l])
        den2 <- (v_val + S_val)^2
        
        b_val <- sum(D[, idx]^2) / (2 * sigma^2)
        w_vec <- c(1, exp(U_safe[t, -l])) / S_val
        term_vec <- D[, -idx] %*% w_vec
        d_val <- sum(term_vec^2) / (2 * sigma^2)
        c_val <- sum(D[, idx] * term_vec) / (2 * sigma^2)
        
        lam1 <- abs(b_val - d_val) * ((v_val^2 * (b_val < d_val)) + (S_val^2 * (b_val > d_val))) / den2
        lam2 <- 2 * (max(b_val, d_val) - c_val) * v_val * S_val / den2
        
        #lam1 <- if(is.finite(lam1)) pmin(max(0, lam1), 500) else 100 
        #lam2 <- if(is.finite(lam2)) pmin(max(0, lam2), 500) else 100
        
        z1 <- rpois(1, lam1); z2 <- rpois(1, lam2)
        b_pg <- 2 * (z1 + z2)
        
        if(b_pg <= 0){
          pseudo_y[t] <- 0; omega[t] <- 0
        } else {
          z_arg <- pmin(pmax(U_safe[t, l] - log(S_val), -20), 20)
          w_val <- max(rpg(1, b_pg, z_arg), 0)
          
          pseudo_y[t] <- log(S_val) + (z1 * (if(b_val < d_val) 1 else -1)) / w_val
          omega[t] <- w_val
        }
      }
      
      # -- FFBS --
      m <- numeric(T+1); C <- numeric(T+1)
      m[1] <- mu; C[1] <- sig_state / max(1 - phi^2, 0) 
      a_hist <- numeric(T); R_hist <- numeric(T)
      
      for(t in 1:T){
        a_t <- phi * m[t] + ar_drift; R_t <- phi^2 * C[t] + sig_state
        a_hist[t] <- a_t; R_hist[t] <- R_t
        
        if(omega[t] < 1e-10){ 
          m[t+1] <- a_t; C[t+1] <- R_t
        } else {
          Q_t <- R_t + 1/omega[t]
          K_t <- R_t / Q_t
          m[t+1] <- a_t + K_t * (pseudo_y[t] - a_t)
          C[t+1] <- R_t * (1 - K_t)
        }
      }
      
      # Backward Sampling
      u_sample <- numeric(T)
      u_curr <- rnorm(1, m[T+1], sqrt(max(C[T+1], 0)))
      u_sample[T] <- u_curr
      for(t in (T-1):1){
        J_t <- C[t+1] * phi / R_hist[t+1]
        h_t <- m[t+1] + J_t * (u_curr - a_hist[t+1])
        H_t <- max(C[t+1] - J_t^2 * R_hist[t+1], 0)
        u_curr <- rnorm(1, h_t, sqrt(H_t)); u_sample[t] <- u_curr
      }
      
      J_0 <- C[1] * phi / R_hist[1]
      h_0 <- m[1] + J_0 * (u_curr - a_hist[1])
      H_0 <- max(C[1] - J_0^2 * R_hist[1], 0)
      new_u0[l] <- rnorm(1, h_0, sqrt(H_0))
      U[, l] <- u_sample
    }
    return(list(U = U, u0 = new_u0))
  }
  
  # --- Initialization & Storage ---
  A <- 0.9 * diag(L-1); B <- rep(0, L-1); u0 <- B; Sigma <- 0.005 * diag(L-1); sigma <- 0.01
  U <- matrix(0, T, L-1)
  
  for (t in 1:T) {
    u1 <- A %*% u0 + t(rmvnorm(1, mean = rep(0, L-1), sigma = Sigma))
    U[t,] <- u1
    u0 <- u1
  }
  u0 <- B
  
  store_U <- array(NA, c(nsim, T, L-1)); store_u0 <- store_A <- store_Sigma <- store_B <- matrix(NA, nsim, L-1); store_sigma <- numeric(nsim)
  
  # --- Main Loop ---
  pb <- txtProgressBar(min = 0, max = burnin + nsim, style = 3)
  for (i in -burnin:nsim) {
    res_AS <- sample_A_Sigma_diag(A, B, Sigma, u0, U)
    A <- res_AS$A; Sigma <- res_AS$Sigma
    B <- sample_B(A, Sigma, u0, U)
    sigma <- sample_sigma(U)
    
    res_U  <- sample_U_scalar_FFBS(U, u0, A, B, Sigma, sigma)
    U <- res_U$U; u0 <- res_U$u0 # Updated jointly
    
    if (i > 0){
      store_U[i,,] <- U; store_u0[i,] <- u0; store_A[i,] <- diag(A)
      store_Sigma[i,] <- diag(Sigma); store_sigma[i] <- sigma; store_B[i,] <- B
    }
    setTxtProgressBar(pb, i + burnin)
  }
  return(list(U = store_U, u0 = store_u0, A = store_A, B = store_B, Sigma = store_Sigma, sigma = store_sigma))
}




# Example Usage
if (FALSE) {
  # Generate dummy data for testing
  T_dummy <- 100
  raw_mat <- matrix(runif(T_dummy * 5), nrow = T_dummy, ncol = 5)
  raw_mat <- t(apply(raw_mat, 1, function(x) x / sum(x))) # Normalize to shares
  cum_data <- t(apply(raw_mat, 1, cumsum))
  Y_dummy <- cum_data[, 1:4] # Use first K columns (K=4, L=5)
  
  # Run MCMC
  res <- run_shape_constrained_fssm(Y_dummy, burnin = 100, nsim = 500)
  
  # Quick Check
  print(colMeans(res$A))
  plot(res$sigma, type='l', main="Traceplot: sigma")
}