library(MCMCpack)
library(statmod)
library(truncnorm)
library(isotone)


###  functional time series with Gaussian process AR(1)  ###
# Y: (T,L)-matrix (T: time series length, L: number of observed points in each time)
# x: L-dimensional vector 
# band: step size of random-walk MH for range parameter in GP
ARGP <- function(Y, x, band=0.05, mc=2000, burn=1000, print=F){
  ## settings
  DD <- as.matrix( dist(x) )
  h.max <- max(DD) 
  TT <- dim(Y)[1]
  L <- dim(Y)[2]
  Ker <- function(mat, h){ exp(-mat/h) }   # kernel function for GP
  
  Inv <- function(Mat){  # matrix inverse for stable computation 
    Mat <- (Mat + t(Mat))/2
    dec <- eigen(Mat)
    lam <- dec$values
    lam[lam<10^(-20)] <- 10^(-20)
    return((dec$vectors)%*%diag(1/lam)%*%t(dec$vectors))
  }
  
  ## objects to store posterior samples 
  Mu.pos <- array(NA, c(mc, TT, L))     # latent functions modeled by GP-AR(1) 
  Sig.pos <- c()   # noise variance 
  h.pos <- c()   # range parameter of GP
  Tau.pos <- c()    # standard deviation of GP 
  rho.pos <- c()   # autocorrelation parameter 
  
  ## prior (for hyperparameter)
  delta <- 0.005   # shape and scale parameters of IG prior for sigma^2 and tau^2
  
  ## initial value
  Mu <- Y
  Sig <- mean((t(Y) - apply(Y, 2, mean))^2)
  h <- median(DD)
  Tau <- sd(apply(Y, 2, mean))
  rho <- 0.5
  
  ## MCMC
  for(k in 1:mc){
    # Mu (latent functions)
    H <- Ker(DD, h)
    IH <- Inv(H)
    for(t in 1:TT){
      # A
      if(t<TT){  A <- Inv( diag(L)/Sig^2 + (1+rho^2)*IH/Tau^2 ) }
      if(t==TT){  A <- Inv( diag(L)/Sig^2 + IH/Tau^2 ) }
      # B
      if(t==1){   B <- as.vector( Y[t,]/Sig^2 + rho*IH%*%Mu[t+1,]/Tau^2 )  }
      if(t>1 & t<TT){   B <- as.vector( Y[t,]/Sig^2 + rho*IH%*%Mu[t+1,]/Tau^2 + rho*IH%*%Mu[t-1,]/Tau^2 )  }
      if(t==TT){   B <- as.vector( Y[t,]/Sig^2 + rho*IH%*%Mu[t-1,]/Tau^2 )  }
      Mu[t,] <- as.vector(mvrnorm(1, A%*%B, A))     
    }
    Mu.pos[k,,] <- Mu
    
    # rho
    A <- 0
    B <- 0
    for(t in 2:TT){
      A <- A + t(Mu[t-1,])%*%IH%*%Mu[t-1,]/Tau^2
      B <- B + t(Mu[t,])%*%IH%*%Mu[t-1,]/Tau^2
    }
    rho <- rtruncnorm(1, a=0, b=1, mean=B/A, sd=sqrt(1/A))
    rho.pos[k] <- rho
    
    # Sigma
    sY <- Y - Mu
    Sig <- sqrt( rinvgamma(1, delta+L*TT/2, delta+sum(sY^2)/2) )
    Sig.pos[k] <- Sig

    # Tau 
    dMu <- matrix(NA, TT, L)
    dMu[1,] <- Mu[1,]
    for(t in 2:TT){
      dMu[t,] <- Mu[t,] - rho*Mu[t-1,]
    }
    mat <- 0
    for(t in 1:TT){
      mat <- mat + t(dMu[t,])%*%IH%*%dMu[t,]
    }
    Tau <- sqrt( rinvgamma(1, delta+L*TT/2, delta+mat/2) )
    Tau.pos[k] <- Tau
    
    # range parameter (h)
    new.h <- h + band*rnorm(1)
    new.h[new.h < 10^(-8)] <- 10^(-8)
    new.H <- Ker(DD, new.h)
    new.IH <- Inv(new.H)
    ss1 <- 0
    ss2 <- 0
    for(t in 1:TT){
      ss1 <- ss1 + t(dMu[t,])%*%IH%*%dMu[t,]/Tau^2
      ss2 <- ss2 + t(dMu[t,])%*%new.IH%*%dMu[t,]/Tau^2
    }
    bb1 <- TT*sum(log(eigen(H)$values))
    bb2 <- TT*sum(log(eigen(new.H)$values))
    val1 <- (-0.5)*bb1 - 0.5*ss1
    val2 <- (-0.5)*bb2 - 0.5*ss2
    prob <- min(1, exp(val2-val1))
    if(runif(1)<prob){
      h <- new.h
    }
    h.pos[k] <- h
    
    # print
    if(print){
      if(round(10*k/mc)==(10*k/mc)){
        print( paste0("MCMC ", round(100*k/mc), "% completed.") )
      }
    } 
  }
  
  # Summary
  om <- 1:burn
  Sig.pos <- Sig.pos[-om]
  Mu.pos <- Mu.pos[-om,,]
  Tau.pos <- Tau.pos[-om]
  h.pos <- h.pos[-om]
  rho.pos <- rho.pos[-om]
  Res <- list(mu=Mu.pos, sig=Sig.pos, tau=Tau.pos, h=h.pos, rho=rho.pos)
  return(Res)
}



###  functional time series with Gaussian process AR(1) with projection  ###
# Y: (T,L)-matrix (T: time series length, L: number of observed points in each time)
# x: L-dimensional vector 
# band: step size of random-walk MH for range parameter in GP
ARGP_proj <- function(Y, x, band=0.05, mc=2000, burn=1000, print=F){
  ## settings
  DD <- as.matrix( dist(x) )
  h.max <- max(DD) 
  TT <- dim(Y)[1]
  L <- dim(Y)[2]
  Ker <- function(mat, h){ exp(-mat/h) }   # kernel function for GP
  
  Inv <- function(Mat){  # matrix inverse for stable computation 
    Mat <- (Mat + t(Mat))/2
    dec <- eigen(Mat)
    lam <- dec$values
    lam[lam<10^(-20)] <- 10^(-20)
    return((dec$vectors)%*%diag(1/lam)%*%t(dec$vectors))
  }
  
  ## objects to store posterior samples 
  Mu.pos <- array(NA, c(mc, TT, L))     # latent functions modeled by GP-AR(1) 
  Sig.pos <- c()   # noise variance 
  h.pos <- c()   # range parameter of GP
  Tau.pos <- c()    # standard deviation of GP 
  rho.pos <- c()   # autocorrelation parameter 
  
  ## prior (for hyperparameter)
  delta <- 0.005   # shape and scale parameters of IG prior for sigma^2 and tau^2
  # 0.02, 0.01, 0.005 
  ## initial value
  Mu <- Y
  Sig <- mean((t(Y) - apply(Y, 2, mean))^2)
  h <- median(DD)
  Tau <- sd(apply(Y, 2, mean))
  rho <- 0.5
  
  ## MCMC
  for(k in 1:mc){
    # Mu (latent functions)
    H <- Ker(DD, h)
    IH <- Inv(H)
    for(t in 1:TT){
      # A
      if(t<TT){  A <- Inv( diag(L)/Sig^2 + (1+rho^2)*IH/Tau^2 ) }
      if(t==TT){  A <- Inv( diag(L)/Sig^2 + IH/Tau^2 ) }
      # B
      if(t==1){   B <- as.vector( Y[t,]/Sig^2 + rho*IH%*%Mu[t+1,]/Tau^2 )  }
      if(t>1 & t<TT){   B <- as.vector( Y[t,]/Sig^2 + rho*IH%*%Mu[t+1,]/Tau^2 + rho*IH%*%Mu[t-1,]/Tau^2 )  }
      if(t==TT){   B <- as.vector( Y[t,]/Sig^2 + rho*IH%*%Mu[t-1,]/Tau^2 )  }
      check <- 1
      while (check >= 1) {
        Mu_candidate <- as.vector(mvrnorm(1, A%*%B, A))
        check <- sum(Mu_candidate < 0)
      }
      Mu[t,] <- ( gpava( c(0,x,1), c(0,Mu_candidate,1) )$x )[c(2:(L+1))]
    }
    Mu.pos[k,,] <- Mu
    
    # rho
    A <- 0
    B <- 0
    for(t in 2:TT){
      A <- A + t(Mu[t-1,])%*%IH%*%Mu[t-1,]/Tau^2
      B <- B + t(Mu[t,])%*%IH%*%Mu[t-1,]/Tau^2
    }
    rho <- rtruncnorm(1, a=0, b=1, mean=B/A, sd=sqrt(1/A))
    rho.pos[k] <- rho
    
    # Sigma
    sY <- Y - Mu
    Sig <- sqrt( rinvgamma(1, delta+L*TT/2, delta+sum(sY^2)/2) )
    Sig.pos[k] <- Sig
    
    # Tau 
    dMu <- matrix(NA, TT, L)
    dMu[1,] <- Mu[1,]
    for(t in 2:TT){
      dMu[t,] <- Mu[t,] - rho*Mu[t-1,]
    }
    mat <- 0
    for(t in 1:TT){
      mat <- mat + t(dMu[t,])%*%IH%*%dMu[t,]
    }
    Tau <- sqrt( rinvgamma(1, delta+L*TT/2, delta+mat/2) )
    Tau.pos[k] <- Tau
    
    # range parameter (h)
    new.h <- h + band*rnorm(1)
    new.h[new.h < 10^(-8)] <- 10^(-8)
    new.H <- Ker(DD, new.h)
    new.IH <- Inv(new.H)
    ss1 <- 0
    ss2 <- 0
    for(t in 1:TT){
      ss1 <- ss1 + t(dMu[t,])%*%IH%*%dMu[t,]/Tau^2
      ss2 <- ss2 + t(dMu[t,])%*%new.IH%*%dMu[t,]/Tau^2
    }
    bb1 <- TT*sum(log(eigen(H)$values))
    bb2 <- TT*sum(log(eigen(new.H)$values))
    val1 <- (-0.5)*bb1 - 0.5*ss1
    val2 <- (-0.5)*bb2 - 0.5*ss2
    prob <- min(1, exp(val2-val1))
    if(runif(1)<prob){
      h <- new.h
    }
    h.pos[k] <- h
    
    # print
    if(print){
      if(round(10*k/mc)==(10*k/mc)){
        print( paste0("MCMC ", round(100*k/mc), "% completed.") )
      }
    } 
  }
  
  # Summary
  om <- 1:burn
  Sig.pos <- Sig.pos[-om]
  Mu.pos <- Mu.pos[-om,,]
  Tau.pos <- Tau.pos[-om]
  h.pos <- h.pos[-om]
  rho.pos <- rho.pos[-om]
  Res <- list(mu=Mu.pos, sig=Sig.pos, tau=Tau.pos, h=h.pos, rho=rho.pos)
  return(Res)
}

