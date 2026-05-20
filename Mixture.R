library(BayesLogit)
library(tidyverse)
library(truncnorm)
library(mvtnorm)
library(MASS)
library(MCMCpack)

run_mixture_fssm <- function(Y, Psi, burnin = 20000, nsim = 80000, seed = 4321) {
  
  set.seed(seed)
  
  # --- Model Setup ---
  T <- nrow(Y)
  K <- ncol(Y)
  L <- ncol(Psi) # Number of basis functions
  
  # Hyperparameters
  a0 <- 0.8; v0 <- 0.2          
  B_a0 <- 0.0; B_v0 <- 5.0      
  und_d <- 0.001; und_s <- 0.001 
  n0 <- 0.001; s0 <- 0.001      
  
  # --- Internal Samplers ---
  
  # 1. Sample Initial State u0
  sample_u0 <- function(U, A, B, Sigma){
    Alpha <- t(t(U) - B)
    Sigma_diag <- diag(Sigma)
    phi_vec <- diag(A)
    
    inv_Sigma0 <- (1 - phi_vec^2) / Sigma_diag
    inv_Sigma  <- 1 / Sigma_diag
    
    prec_post <- inv_Sigma0 + (phi_vec^2) * inv_Sigma
    var_post  <- 1 / prec_post
    mean_post <- var_post * (phi_vec * inv_Sigma * Alpha[1,])
    
    alpha0 <- rnorm(L-1, mean = mean_post, sd = sqrt(var_post))
    return(alpha0 + B)
  }
  
  # 2. Sample AR(1) coef A and State Variance Sigma
  sample_A_Sigma_diag <- function(A, B, Sigma, u0, U){
    Alpha  <- t(t(U) - B)
    alpha0 <- u0 - B
    
    new_A     <- numeric(L-1)
    new_Sigma <- numeric(L-1)
    
    x0_mat <- rbind(alpha0, Alpha[1:(T-1), , drop=FALSE])
    x1_mat <- Alpha
    
    for (l in 1:(L-1)){
      x0 <- x0_mat[,l]
      x1 <- x1_mat[,l]
      sig_sq <- Sigma[l,l]
      
      # Sample A (Metropolis step)
      denom    <- sum(x0^2) / sig_sq + 1/v0^2
      var_phi  <- 1 / denom
      mean_phi <- var_phi * (sum(x0 * x1) / sig_sq + a0/v0^2)
      
      old_phi  <- A[l,l]
      cand_phi <- rtruncnorm(1, a=-1, b=1, mean=mean_phi, sd=sqrt(var_phi))
      
      ratio <- sqrt((1 - cand_phi^2) / (1 - old_phi^2))
      if(runif(1) < ratio){
        new_A[l] <- cand_phi
      } else {
        new_A[l] <- old_phi
      }
      
      # Sample Sigma
      phi_curr <- new_A[l]
      resid <- x1 - phi_curr * x0
      shape_post <- (T + und_d + 1) / 2
      rate_post  <- (sum(resid^2) + und_s + (1 - phi_curr^2) * alpha0[l]^2) / 2
      new_Sigma[l] <- 1 / rgamma(1, shape = shape_post, rate = rate_post)
    }
    
    return(list(A = diag(new_A), Sigma = diag(new_Sigma)))
  }
  
  # 3. Sample Mean Level B
  sample_B <- function(A, Sigma, u0, U){
    sampleB <- numeric(L-1)
    U_lag <- rbind(u0, U[1:(T-1), , drop=FALSE])
    
    for (l in 1:(L-1)){
      phi <- A[l,l]
      sig <- Sigma[l,l]
      
      y_vec <- U[,l] - phi * U_lag[,l]
      x_val <- (1 - phi)
      
      prec_prior <- 1 / B_v0^2
      prec_like  <- T * x_val^2 / sig
      
      var_post  <- 1 / (prec_prior + prec_like)
      mean_post <- var_post * ( (sum(y_vec) * x_val / sig) + (B_a0 / B_v0^2) )
      
      # u0 correction
      prec_u0 <- (1 - phi^2) / sig
      var_exact <- 1 / (1/var_post + prec_u0)
      mean_exact <- var_exact * (mean_post/var_post + u0[l] * prec_u0)
      
      sampleB[l] <- rnorm(1, mean_exact, sqrt(var_exact))
    }
    return(sampleB)
  }
  
  # 4. Sample Latent Indicators Z (修正: sigmas は長さ L のベクトル)
  sample_Z <- function(U, sigmas){
    result_Z <- matrix(0, nrow = T, ncol = K)
    exp_U <- exp(pmin(pmax(U, -20), 20))
    Pi <- cbind(1, exp_U) / (1 + rowSums(exp_U)) 
    
    for(t in 1:T){
      for(k in 1:K){
        log_lik <- -log(sigmas) - 0.5 * (Y[t,k] - Psi[k,])^2 / (sigmas^2)
        log_prob <- log(Pi[t,]) + log_lik
        log_prob <- log_prob - max(log_prob)
        prob <- exp(log_prob)
        result_Z[t,k] <- sample(1:L, size=1, prob=prob/sum(prob))
      }
    }
    return(result_Z)
  }
  
  # 5. Sample Observation Variances
  sample_sigmas_conj <- function(Z){
    new_sigmas <- numeric(L)
    for(l in 1:L){
      indices <- which(Z == l, arr.ind = TRUE)
      n_l <- nrow(indices)
      
      ssr <- 0
      if(n_l > 0){
        for(idx in 1:n_l){
          t_idx <- indices[idx, 1]
          k_idx <- indices[idx, 2]
          ssr <- ssr + (Y[t_idx, k_idx] - Psi[k_idx, l])^2
        }
      }
      shape_post <- (2.0 + n_l) / 2
      rate_post  <- (0.01 + ssr) / 2
      new_sigmas[l] <- sqrt(1 / rgamma(1, shape = shape_post, rate = rate_post))
    }
    return(new_sigmas)
  }
  
  # 6. Sample States U
  sample_U_mix_FFBS <- function(U, u0, Z, A, B, Sigma){
    
    for (l in 1:(L-1)) {
      phi <- A[l, l]
      sig_state <- Sigma[l, l]
      mu  <- B[l]
      ar_drift <- mu * (1 - phi)
      
      pseudo_y <- numeric(T)
      omega    <- numeric(T)
      
      # Sum of exp of other states
      if ((L-1) == 1) {
        S_vec <- rep(1, T)
      } else {
        S_vec <- 1 + rowSums(exp(U[, -l, drop=FALSE]))
      }
      
      # -- PG Augmentation --
      for (t in 1:T) {
        # Z maps to 1..L. Current U column corresponds to class l+1
        target_class <- l + 1
        
        n_success <- sum(Z[t,] == target_class)
        n_trials  <- K 

        psi_val <- U[t, l] - log(S_vec[t])
        
        # Sample PG
        w_val <- sum(rpg(num = K, h = 1, z = psi_val))
        w_val <- max(w_val, 1e-9)
        
        kappa <- n_success - n_trials / 2
        
        pseudo_y[t] <- kappa / w_val + log(S_vec[t])
        omega[t]    <- w_val
      }
      
      # -- FFBS --
      m <- numeric(T+1); C <- numeric(T+1)
      m[1] <- u0[l]; C[1] <- 0
      
      a_hist <- numeric(T); R_hist <- numeric(T)
      
      # Forward Filter
      for(t in 1:T){
        a_t <- phi * m[t] + ar_drift
        R_t <- phi^2 * C[t] + sig_state
        
        a_hist[t] <- a_t; R_hist[t] <- R_t
        
        Q_t <- R_t + 1/omega[t]
        K_t <- R_t / Q_t
        
        m[t+1] <- a_t + K_t * (pseudo_y[t] - a_t)
        C[t+1] <- R_t * (1 - K_t)
      }
      
      # Backward Sampler
      u_sample <- numeric(T)
      u_next   <- rnorm(1, m[T+1], sqrt(C[T+1]))
      u_sample[T] <- u_next
      
      for(t in (T-1):1){
        J_t <- C[t+1] * phi / R_hist[t+1]
        h_t <- m[t+1] + J_t * (u_next - a_hist[t+1])
        H_t <- C[t+1] - J_t^2 * R_hist[t+1]
        
        if(H_t < 0) H_t <- 1e-9
        
        u_next <- rnorm(1, h_t, sqrt(H_t))
        u_sample[t] <- u_next
      }
      
      U[, l] <- u_sample
    }
    return(list(U = U))
  }
  
  # --- Initialization ---
  A     <- 0.9 * diag(L-1)
  B     <- rep(0, L-1)
  u0    <- B
  Sigma <- 0.005 * diag(L-1)
  
  # Initialize Mixture variances
  sigmas <- rep(0.01, L)
  
  # Initialize U
  U <- matrix(0, nrow = T, ncol = L-1)
  u_curr <- u0
  for (t in 1:T) {
    U[t,] <- A %*% u_curr + rnorm(L-1, 0, sqrt(diag(Sigma)))
    u_curr <- U[t,]
  }
  u0 <- B
  
  # Initialize Z 
  Z <- matrix(sample(1:L, T*K, replace=TRUE), nrow=T, ncol=K)
  
  # Storage
  store_U      <- array(NA, dim = c(nsim, T, L-1))
  store_u0     <- matrix(NA, nrow = nsim, ncol = L-1)
  store_A      <- matrix(NA, nrow = nsim, ncol = L-1)
  store_Sigma  <- matrix(NA, nrow = nsim, ncol = L-1)
  store_sigmas <- matrix(NA, nrow = nsim, ncol = L)
  store_B      <- matrix(NA, nrow = nsim, ncol = L-1)
  
  # --- Main MCMC Loop ---
  cat("Starting MCMC (Mixture)...\n")
  pb <- txtProgressBar(min = 0, max = burnin + nsim, style = 3)
  
  for (i in -burnin:nsim) {
    
    # 1. Update Hyperparameters
    res_AS <- sample_A_Sigma_diag(A, B, Sigma, u0, U)
    A      <- res_AS$A
    Sigma  <- res_AS$Sigma
    B      <- sample_B(A, Sigma, u0, U)
    u0     <- sample_u0(U, A, B, Sigma)
    
    # 2. Update Latent Classes Z
    Z <- sample_Z(U, sigmas)
    
    # 3. Update Variances sigmas
    sigmas <- sample_sigmas_conj(Z)
    
    # 4. Update States U
    res_U  <- sample_U_mix_FFBS(U, u0, Z, A, B, Sigma)
    U      <- res_U$U
    
    # 5. Store
    if (i > 0){
      store_U[i,,] <- U
      store_u0[i,] <- u0
      store_A[i,]  <- diag(A)
      store_Sigma[i,] <- diag(Sigma)
      store_sigmas[i,] <- sigmas
      store_B[i,]     <- B
    }
    
    setTxtProgressBar(pb, i + burnin)
  }
  close(pb)
  
  return(list(U = store_U, u0 = store_u0, A = store_A, B = store_B, 
              Sigma = store_Sigma, sigmas = store_sigmas))
}






# Example Usage

if (FALSE) {
  # Generate dummy data
  T_dummy <- 100
  raw_mat <- matrix(runif(T_dummy * 5), nrow = T_dummy, ncol = 5)
  raw_mat <- t(apply(raw_mat, 1, function(x) x / sum(x))) 
  cum_data <- t(apply(raw_mat, 1, cumsum))
  Y_dummy <- cum_data[, 1:4] # Use first K columns
  
  # Run MCMC
  res <- run_mixture_fssm(Y_dummy, burnin = 100, nsim = 200)
  
  # Quick Check
  print(colMeans(res$A))
  # Check convergence of a specific variance parameter
  plot(res$sigmas[,1,1], type='l', main="Traceplot: sigma[1,1]")
}