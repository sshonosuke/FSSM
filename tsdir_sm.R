library(MCMCpack)
library(mvtnorm)
library(truncnorm)
library(MASS)

run_sm_tsdir <- function(Y, burnin = 1000, nsim = 2000, seed = 4321, 
                         rw_step_h = 0.1, rw_step_lam = 0.5) {
  
  set.seed(seed)
  
  # Data Preparation
  T_len <- nrow(Y)
  K <- ncol(Y) 
  
  # Convert cumulative Y to interval proportions
  Y_aug <- cbind(0, Y, 1)
  Y_prop <- t(apply(Y_aug, 1, diff))
  
  # Define grid points x 
  x_points <- seq(1/(K+1), 1-1/(1+K), length.out = K)
  
  ### Functions ###
  
  LorenzSM <- function(p_vec, a, q) {
    if (a <= 0 || q <= 0 || (a*q) <= 1) return(rep(NA, length(p_vec)))
    shape1 <- 1 + 1/a
    shape2 <- q - 1/a
    z <- 1 - (1 - p_vec)^(1/q)
    return(pbeta(z, shape1, shape2))
  }
  
  GiniSM <- function(a, q) {
    if (a <= 0 || q <= 0 || (a*q) <= 1) return(NA)
    val <- lgamma(q) + lgamma(2*q - 1/a) - lgamma(q - 1/a) - lgamma(2*q)
    return(1 - exp(val))
  }
  
  log_lik_dirichlet <- function(y_row, theta, lambda) {
    a <- theta[1]
    q <- theta[2]
    l_vals <- LorenzSM(x_points, a, q)
    if (any(is.na(l_vals))) return(-Inf)
    l_aug <- c(0, l_vals, 1)
    l_diff <- diff(l_aug)
    alpha <- lambda * l_diff
    if (any(alpha <= 0)) return(-Inf)
    ll <- lgamma(sum(alpha)) - sum(lgamma(alpha)) + sum((alpha - 1) * log(y_row))
    return(ll)
  }
  
  # MCMC
  D <- 2
  H <- matrix(0, nrow = T_len, ncol = D)
  H[,1] <- log(3.0) 
  H[,2] <- log(1.5) 
  
  # Parameters
  mu <- colMeans(H)
  Phi <- diag(0.9, D) 
  Sigma <- diag(0.01, D)
  lambda <- 100 
  
  # Priors
  mu0 <- c(0, 0); Sigma0_mu <- diag(10, D)
  S0_wish <- diag(0.1, D); v0_wish <- D + 2
  
  # Storage
  store_H <- array(NA, dim = c(nsim, T_len, D))
  store_mu <- matrix(NA, nrow = nsim, ncol = D)
  store_Phi <- matrix(NA, nrow = nsim, ncol = D)
  store_Sigma <- array(NA, dim = c(nsim, D, D))
  store_lambda <- numeric(nsim)
  store_Gini <- matrix(NA, nrow = nsim, ncol = T_len)
  
  # Main Loop
  cat("Starting MCMC for Singh-Maddala Model...\n")
  pb <- txtProgressBar(min = 0, max = burnin + nsim, style = 3)
  
  for (iter in -burnin:nsim) {
    
    # Pre-compute inverse variance for AR part
    Sigma_inv <- solve(Sigma)
    
    # Sample H_t
    for (t in 1:T_len) {
      h_curr <- H[t, ]
      h_prop <- h_curr + rnorm(D, 0, rw_step_h)
      theta_curr <- exp(h_curr)
      theta_prop <- exp(h_prop)
      
      if (theta_prop[1] * theta_prop[2] > 1) {
        
        ll_curr <- log_lik_dirichlet(Y_prop[t, ], theta_curr, lambda)
        ll_prop <- log_lik_dirichlet(Y_prop[t, ], theta_prop, lambda)
        
        # calculate AR log-prior
        calc_lp_ar <- function(val, t_idx) {
          lp <- 0
          if (t_idx > 1) {
            mean_prev <- mu + as.vector(Phi %*% (H[t_idx-1, ] - mu))
            res <- val - mean_prev
            lp <- lp - 0.5 * t(res) %*% Sigma_inv %*% res
          } else {
            lp <- lp - 0.5 * sum((val - mu)^2) 
          }
          if (t_idx < T_len) {
            mean_next <- mu + as.vector(Phi %*% (val - mu))
            res_next <- H[t_idx+1, ] - mean_next
            lp <- lp - 0.5 * t(res_next) %*% Sigma_inv %*% res_next
          }
          return(as.numeric(lp))
        }
        
        lp_ar_curr <- calc_lp_ar(h_curr, t)
        lp_ar_prop <- calc_lp_ar(h_prop, t)
        
        log_ratio <- (ll_prop + lp_ar_prop) - (ll_curr + lp_ar_curr)
        
        if (log(runif(1)) < log_ratio) {
          H[t, ] <- h_prop
        }
      }
    }
    
    # Sample Lambda
    lam_prop <- exp(log(lambda) + rnorm(1, 0, rw_step_lam))
    ll_lam_curr <- 0; ll_lam_prop <- 0
    valid_lam <- TRUE
    
    for(t in 1:T_len) {
      lc <- log_lik_dirichlet(Y_prop[t,], exp(H[t,]), lambda)
      lp <- log_lik_dirichlet(Y_prop[t,], exp(H[t,]), lam_prop)
      if(lp == -Inf) valid_lam <- FALSE
      ll_lam_curr <- ll_lam_curr + lc
      ll_lam_prop <- ll_lam_prop + lp
    }
    
    if (valid_lam) {
      if (log(runif(1)) < (ll_lam_prop - ll_lam_curr)) {
        lambda <- lam_prop
      }
    }
    
    # Sample Mu, Phi, Sigma
    
    # 1 Sample Mu
    mu_hat <- colMeans(H)
    mu <- as.vector(rmvnorm(1, mu_hat, Sigma / T_len))
    
    # 2 Sample Phi & Sigma
    X_reg <- H[1:(T_len-1), ] - matrix(rep(mu, T_len-1), ncol=D, byrow=TRUE)
    Y_reg <- H[2:T_len, ] - matrix(rep(mu, T_len-1), ncol=D, byrow=TRUE)
    
    S_post <- S0_wish + t(Y_reg - Y_reg %*% diag(diag(Phi))) %*% (Y_reg - Y_reg %*% diag(diag(Phi)))
    Sigma <- riwish(v0_wish + T_len, S_post)
    
    for(d in 1:D) {
      num <- sum(X_reg[,d] * Y_reg[,d])
      den <- sum(X_reg[,d]^2)
      phi_mean <- num / den
      phi_var <- Sigma[d,d] / den
      
      phi_cand <- rtruncnorm(1, a=-0.99, b=0.99, mean=phi_mean, sd=sqrt(phi_var))
      Phi[d,d] <- phi_cand
    }
    
    # Store
    if (iter > 0) {
      store_H[iter, , ] <- H
      store_mu[iter, ] <- mu
      store_Phi[iter, ] <- diag(Phi)
      store_Sigma[iter, , ] <- Sigma
      store_lambda[iter] <- lambda
      
      for(t in 1:T_len) {
        theta <- exp(H[t,])
        store_Gini[iter, t] <- GiniSM(theta[1], theta[2])
      }
    }
    setTxtProgressBar(pb, iter + burnin)
  }
  close(pb)
  
  return(list(
    H = store_H, Theta = exp(store_H), mu = store_mu,
    Phi = store_Phi, Sigma = store_Sigma, lambda = store_lambda,
    Gini = store_Gini
  ))
}







# --- Example Usage ---
if(FALSE){
  res_sm <- run_sm_tsdir(Y, burnin=1000, nsim=2000)
  
  # Plot Gini
  gini_mean <- colMeans(res_sm$Gini, na.rm=TRUE)
  plot(gini_mean, type='l', main="Estimated Gini (SM Model)")
}
