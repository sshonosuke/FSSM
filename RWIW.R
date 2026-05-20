library(dlm)
library(MASS)
library(isotone)

run_rwiw_dlm <- function(Y, burnin = 20000, nsim = 80000, seed = 4321) {
  
  set.seed(seed)
  
  # --- Model Setup ---
  TT <- nrow(Y)
  m  <- ncol(Y) 
  x_points <- seq(1/(1+m), 1-1/(1+m), length.out = m)
  
  # Hyperparameters
  delta0 <- 1; delta1 <- 1
  V0 <- (delta0 + 4 - 1) * diag(0.02^2, m)
  W0 <- (delta1 + 4 - 1) * diag(0.01^2, m)
  
  # Define DLM structure
  mod <- dlm(FF = diag(m),
             V  = diag(m),
             GG = diag(m),
             W  = diag(m),
             m0 = rep(0, m),
             C0 = diag(100^2, m))
  
  # Initial values
  mod$V <- V0 / (delta0 + m - 1)
  mod$W <- W0 / (delta1 + m - 1)
  
  # Storage
  store_Theta <- array(NA, dim=c(nsim, TT+1, m)) 
  store_V     <- array(NA, dim=c(nsim, m, m))
  store_W     <- array(NA, dim=c(nsim, m, m))
  
  # --- MCMC Loop ---
  cat("Starting Gibbs Sampling (RWIW)...\n")
  pb <- txtProgressBar(min = 0, max = burnin + nsim, style = 3)
  
  for(i in -burnin:nsim){
    
    # 1. Sample States (FFBS)
    modFilt <- dlmFilter(Y, mod, simplify=TRUE)
    
    # Positivity check
    check <- 1
    count_check <- 0
    theta_raw <- NULL
    
    while (check >= 1) {
      theta_raw <- dlmBSample(modFilt)
      check <- sum(theta_raw < 0)
      count_check <- count_check + 1
      if(count_check > 50) {
        check <- 0 
      }
    }
    
    # PAVA
    theta_curr <- matrix(0, nrow=TT+1, ncol=m)
    for(t in 1:(TT+1)){
      pava_res <- gpava(c(0, x_points, 1), c(0, theta_raw[t,], 1))
      theta_curr[t,] <- pava_res$x[2:(m+1)]
    }
    
    # Sample V
    res_y <- Y - theta_curr[-1,] 
    S_post <- crossprod(res_y) + V0
    
    prec_V <- rwishart(df = delta0 + 1 + TT, p = m, Sigma = solve(S_post))
    mod$V <- solve(prec_V)
    
    # Sample W
    theta_center <- theta_curr[-1,] - theta_curr[-(TT+1),]
    SS_post <- crossprod(theta_center) + W0
    
    prec_W <- rwishart(df = delta1 + 1 + TT, p = m, Sigma = solve(SS_post))
    mod$W <- solve(prec_W)
    
    # Store Samples
    if(i > 0){
      store_Theta[i,,] <- theta_curr
      store_V[i,,]     <- mod$V
      store_W[i,,]     <- mod$W
    }
    
    setTxtProgressBar(pb, i + burnin)
  }
  close(pb)
  
  return(list(
    Theta = store_Theta,
    V     = store_V,
    W     = store_W
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
  
  # Run MCMC
  res <- run_rwiw_dlm(Y, burnin = 100, nsim = 200)
  
  # Check results
  print(dim(res$Theta))
}