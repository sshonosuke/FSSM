library(BayesLogit); library(tidyverse); library(truncnorm); library(mvtnorm)
library(MASS); library(MCMCpack); library(isotone); library(coda)

source("FSSM.R"); source("Mixture.R")
source("RWIW.R"); source("RWSV.R"); source("tsdir_sm.R")
source("ARGP-function.R")

TT <- 200; burnin <- 5000; nsim <- 20000
x_all <- seq(0.05, 0.95, by = 0.05)
M_fine <- 100
x_fine <- seq(0, 1, length.out = M_fine + 1)

# Rasche型
true_ftx <- function(x, t, T_max) {
  a <- 0.83 + 0.01*sin(2*pi*t/T_max)
  b  <- 0.73 + 0.01*cos(2*pi*t/T_max)
  return( (1-(1-x)^a)^{1/b} )
}
true_Gini <- function(t, T_max) {
  a <- 0.83 + 0.01*sin(2*pi*t/T_max)
  b  <- 0.73 + 0.01*cos(2*pi*t/T_max)
  return(1 - 2*beta(1/a, 1/b+1)/a)
}
Gini_beta <- function(a, b) { 1 - 2 * b / (a + b) }
calc_gini_trapezoidal <- function(x_pts, y_pts) {
  x_a <- c(0, x_pts, 1); y_a <- c(0, y_pts, 1)
  area <- sum(diff(x_a) * (head(y_a, -1) + tail(y_a, -1)) / 2)
  return(1 - 2 * area)
}
LorenzSM_pbeta <- function(p, a, q) {
  z <- 1-(1-p)^(1/q); pbeta(z, 1+1/a, q-1/a)
}

# basis set
basis_sets <- list(
  Basis1 = list(a = c(1.0, 1.5, 3.0, 1.0, 1.0), b = c(1.0, 1.0, 1.0, 0.7, 0.3)), 
  Basis3 = list(a = c(1.2, 1.5, 1.0, 1.3, 1.3), b = c(0.9, 0.8, 0.6, 0.8, 0.7))                       # L=3
)

visual_check_basis <- function(K_val, Y_obs, x_obs, basis_sets, TT) {
  par(mfrow = c(1, length(basis_sets)), mar = c(4, 4, 3, 1))
  x_f <- seq(0, 1, length.out = 100)
  col_t <- rgb(0.5, 0, 0, 0.05)
  col_o <- rgb(0, 0, 1, 0.1)
  
  for (b_n in names(basis_sets)) {
    bs <- basis_sets[[b_n]]
    plot(NULL, xlim=c(0,1), ylim=c(0,1), 
         main=paste0("K=", K_val, " / ", b_n),
         xlab="x", ylab="f(x)")
    # basis plot
    for (l in 1:length(bs$a)) lines(x_f, pbeta(x_f, bs$a[l], bs$b[l]), col="skyblue")
    
    for (t in 1:TT) {
      lines(x_f, true_ftx(x_f, t, TT), col=col_t)
      points(x_obs, Y_obs[t,], col=col_o, pch=16, cex=0.5)
    }
  }
}

predict_gp <- function(x_new, x_obs, y_obs, h, tau, rho, sig) {
  
  calc_kernel <- function(x1, x2, h, tau) {
    d_mat <- abs(outer(x1, x2, "-"))
    return(tau^2 * exp(-d_mat / h))
  }
  
  K_xx <- calc_kernel(x_obs, x_obs, h, tau)
  
  diag(K_xx) <- diag(K_xx) + max(sig^2, 1e-9)
  
  K_sx <- calc_kernel(x_new, x_obs, h, tau)
  
  f_mean <- K_sx %*% solve(K_xx, y_obs)
  
  return(as.numeric(f_mean))
}

true_g_series <- sapply(1:TT, function(t) true_Gini(t, TT))

plot(true_g_series, type = "l", lwd = 2, col = "black",
     main = "True Gini Coefficient Time Series",
     xlab = "Time (t)", ylab = "Gini Index",
     ylim = c(min(true_g_series)-0.05, max(true_g_series)+0.05))

par(mfrow = c(2, 1), mar = c(4, 4, 2, 8), xpd = TRUE)
res_storage <- list() # サマリー用

for (K in c(9)) {
  # DGP
  set.seed(42) 
  x_obs <- seq(1/(K+1), K/(K+1), length.out = K)
  Y_obs <- matrix(0, TT, K)
  true_g <- numeric(TT)
  for (t in 1:TT) {
    Y_obs[t,] <- true_ftx(x_obs, t, TT) + rnorm(K, 0, 0.005)* sqrt(1 - 4*(x_obs-0.5)^2)
    Y_obs[t,] <- cumsum(abs(diff(c(0, Y_obs[t,]))))[1:K]
    true_g[t] <- true_Gini(t, TT)
  }
  
  visual_check_basis(K, Y_obs, x_obs, basis_sets, TT)
  
  cat(sprintf("\nProcessing K = %d...\n", K))
  
  # FSSM (L=5, Basis1)
  Psi_o1 <- matrix(0, K, 5); for(l in 1:5) Psi_o1[,l] <- pbeta(x_obs, basis_sets$Basis1$a[l], basis_sets$Basis1$b[l])
  res_f1 <- run_shape_constrained_fssm(Y_obs, Psi_o1, burnin=burnin, nsim=nsim)
  g_f1 <- rowMeans(matrix(sapply(1:nsim, function(i) {
    pi_t <- cbind(1, exp(res_f1$U[i,,])) / (1 + rowSums(exp(res_f1$U[i,,])))
    pi_t %*% Gini_beta(basis_sets$Basis1$a, basis_sets$Basis1$b)
  }), nrow=TT))
  if(K==9) res_storage$f1 <- res_f1
  
  # FSSM (L=5, Basis3)
  Psi_o3 <- matrix(0, K, 5); for(l in 1:5) Psi_o3[,l] <- pbeta(x_obs, basis_sets$Basis3$a[l], basis_sets$Basis3$b[l])
  res_f3 <- run_shape_constrained_fssm(Y_obs, Psi_o3, burnin=burnin, nsim=nsim)
  g_f3 <- rowMeans(matrix(sapply(1:nsim, function(i) {
    pi_t <- cbind(1, exp(res_f3$U[i,,])) / (1 + rowSums(exp(res_f3$U[i,,])))
    pi_t %*% Gini_beta(basis_sets$Basis3$a, basis_sets$Basis3$b)
  }), nrow=TT))
  
  # ARGP (Projected)
  res_a <- ARGP_proj(Y_obs, x_obs, mc=nsim+burnin, burn=burnin)
  g_a <- numeric(TT)
  for(t in 1:TT) {
    f_fine <- predict_gp(x_fine, c(0, x_obs, 1), c(0, colMeans(res_a$mu[,t,]), 1), mean(res_a$h), mean(res_a$tau), 0, mean(res_a$sig))
    g_a[t] <- calc_gini_trapezoidal(x_fine[2:M_fine], f_fine[2:M_fine])
  }
  res_a <- 0
  
  # DLM (RWSV)
  res_sv <- run_rwsv_dlm(Y_obs, nsim=nsim) 
  g_sv <- res_sv$gini
  res_sv <- 0
  
  # TS-DIR
  res_k <- run_sm_tsdir(Y_obs, burnin=burnin, nsim=nsim)
  g_k <- colMeans(res_k$Gini)
  res_k <- 0
  
  # results
  # plots
  par(mfrow = c(2, 1), mar = c(4, 4, 3, 10), xpd = TRUE)
  
  plot_cols <- c(rgb(0, 0, 0, alpha = 0.3), "blue", "red", "darkgreen", "purple", "orange")
  plot_ltys <- c(NA, 1, 2, 3, 4, 5)

  plot(true_g, type="p", pch=20, col=plot_cols[1], cex=0.6, 
       ylim=c(min(true_g)-0.02, max(true_g)+0.02),
       main=paste("Time-series of Gini Coefficients (K =", K, ")"), 
       xlab="Time", ylab="Gini Index")
  
  lines(g_f1, col=plot_cols[2], lwd=1.5, lty=plot_ltys[2]) # Basis 1
  lines(g_f3, col=plot_cols[3], lwd=1.5, lty=plot_ltys[3]) # Basis 3
  lines(g_a,  col=plot_cols[4], lwd=1.5, lty=plot_ltys[4]) # ARGP
  lines(g_sv, col=plot_cols[5], lwd=1.5, lty=plot_ltys[5]) # DLM(SV)
  lines(g_k,  col=plot_cols[6], lwd=1.5, lty=plot_ltys[6]) # TS-DIR
  
  # legend
  legend("topright", inset=c(-0.42, 0), bty="n", 
         legend=c("True", "Basis 1", "Basis 3", "ARGP(projected)", "DLM(SV)", "TS-DIR"),
         col=plot_cols, lty=plot_ltys, pch=c(20, rep(NA, 5)), lwd=1.5, cex=0.8)
  
  write.csv(cbind(true_g, g_f1, g_f3, g_a, g_sv, g_k), file=paste0("ginis_K_", K, ".csv"))
}








summarize_fssm_corrected <- function(res) {
  means <- c(apply(res$A, 2, mean), apply(res$Sigma, 2, mean), apply(res$B, 2, mean), mean(res$sigma))
  sds   <- c(apply(res$A, 2, sd),   apply(res$Sigma, 2, sd),   apply(res$B, 2, sd),   sd(res$sigma))
  
  # 95% CI
  quants <- rbind(
    t(apply(res$A, 2, quantile, c(0.025, 0.975))),
    t(apply(res$Sigma, 2, quantile, c(0.025, 0.975))),
    t(apply(res$B, 2, quantile, c(0.025, 0.975))),
    quantile(res$sigma, c(0.025, 0.975))
  )
  
  # ESS
  ess_vals <- c(
    apply(res$A, 2, function(x) effectiveSize(as.mcmc(x))),
    apply(res$Sigma, 2, function(x) effectiveSize(as.mcmc(x))),
    apply(res$B, 2, function(x) effectiveSize(as.mcmc(x))),
    effectiveSize(as.mcmc(res$sigma))
  )
  
  summary_df <- data.frame(
    Mean = means,
    SD   = sds,
    `Q2.5` = quants[,1],
    `Q97.5` = quants[,2],
    ESS  = round(ess_vals)
  )
  
  rownames(summary_df) <- c("phi_1", "phi_2", "phi_3", "phi_4", "v2_1", "v2_2", "v2_3", "v2_4", "mu_1", "mu_2", "mu_3", "mu_4", "sigma_obs")
  
  return(round(summary_df, 4))
}

cat("\n--- Table: Corrected Posterior Summary (FSSM L=5, K=4) ---\n")
print(summarize_fssm_corrected(res_storage$f1))


#res_test <- res_storage$f1
#plot(res_test$B[,4], type="l")


# pi plots
plot_pi_evolution <- function(res, L, title, colors) {

  pi_mean <- matrix(0, TT, L)
  for(i in 1:nsim) {
    u_i <- res$U[i,,]
    exp_u <- exp(u_i)
    denom <- 1 + rowSums(exp_u)
    pi_mean <- pi_mean + cbind(1/denom, exp_u/denom)
  }
  pi_mean <- pi_mean / nsim
  
  pi_cum <- t(apply(pi_mean, 1, cumsum))
  
  plot(NULL, xlim=c(1, TT), ylim=c(0, 1), main=title, 
       xlab="Time", ylab="Weights (pi)", xaxs="i", yaxs="i")
  
  for (l in L:1) {
    upper <- pi_cum[, l]
    lower <- if(l == 1) rep(0, TT) else pi_cum[, l-1]
    polygon(c(1:TT, TT:1), c(upper, rev(lower)), 
            col = colors[l], border = "white", lwd = 0.5)
  }
  
  # legend
  legend("right", inset=c(-0.25, 0), bty="n", 
         legend=paste0("Basis ", 1:L), fill=colors, cex=0.7)
}


par(mfrow = c(3, 1), mar = c(4, 4, 3, 10), xpd = TRUE)

# Basis 1 (L=5)
plot_pi_evolution(res_f1, L=5, title="Weight Evolution (Basis 1)", 
                  colors=rev(terrain.colors(5)))

# Basis 3 (L=5)
plot_pi_evolution(res_f3, L=5, title="Weight Evolution (Basis 3)", 
                  colors=rev(terrain.colors(5)))

#plot( colMeans(res_f1$U[,,3]), type="l" )
#effectiveSize( as.mcmc(res_f1$U[,,4]) )







# 1. Load Data
df4 <- read.csv("ginis_K_4.csv")
df9 <- read.csv("ginis_K_9.csv")

# 2. Plot Settings
#plot_cols <- c(rgb(0, 0, 0, alpha = 0.3), "blue", "red", "darkgreen", "purple", "orange")
plot_cols <- c(rgb(0, 0, 0, alpha = 0.5), "red", "blue", "darkgreen", "purple", "orange")
plot_ltys <- c(NA, 1, 2, 4, 5, 6) 
plot_labels <- c("True", "Basis 1", "Basis 3", "ARGP", "DLM (SV)", "TS-DIR")

# 3. Consistent Y-axis range
all_values <- unlist(df4[, -1], df9[, -1])
ylim_range <- range(all_values, na.rm = TRUE)
ylim_range[1] <- ylim_range[1] - 0.003
ylim_range[2] <- ylim_range[2] + 0.003

# 4. Save Settings (Reduced height to 400 for a tighter look)
#png("Figure1_Gini_Comparison.png", width = 1000, height = 400, res = 120)
pdf("Figure1_Gini_Comparison.pdf", width = 6, height = 3)

# 5. Layout Setup
layout_mat <- matrix(c(1, 2, 3, 3), nrow = 2, ncol = 2, byrow = TRUE)
layout(layout_mat, heights = c(10, 1.2))

# oma: Reduced top margin (3rd value) to 0.5
par(oma = c(0, 1, 0.5, 1))

# --- Left Plot (K = 4) ---
# mar: Reduced top margin (3rd value) to 1.2
par(mar = c(2, 4.8, 1.8, 0.2)) 
plot(df4$true_g, type = "p", pch = 20, col = plot_cols[1], cex = 0.3,
     ylim = ylim_range, xlab = "", ylab = "", main = "K = 4", 
     yaxt = "s")

mtext("Gini Coefficient", side = 2, line = 3.2, outer = FALSE, cex = 0.9)

lines(df4$g_f1, col = plot_cols[2], lwd = 1.2, lty = plot_ltys[2])
lines(df4$g_f3, col = plot_cols[3], lwd = 1.2, lty = plot_ltys[3])
lines(df4$g_a,  col = plot_cols[4], lwd = 1.2, lty = plot_ltys[4])
lines(df4$g_sv, col = plot_cols[5], lwd = 1.2, lty = plot_ltys[5])
lines(df4$g_k,  col = plot_cols[6], lwd = 1.2, lty = plot_ltys[6])

# --- Right Plot (K = 9) ---
# Symmetric margins
par(mar = c(2, 0.2, 1.8, 4.8))
plot(df9$true_g, type = "p", pch = 20, col = plot_cols[1], cex = 0.3,
     ylim = ylim_range, xlab = "", ylab = "", main = "K = 9", 
     yaxt = "n")

lines(df9$g_f1, col = plot_cols[2], lwd = 1.2, lty = plot_ltys[2])
lines(df9$g_f3, col = plot_cols[3], lwd = 1.2, lty = plot_ltys[3])
lines(df9$g_a,  col = plot_cols[4], lwd = 1.2, lty = plot_ltys[4])
lines(df9$g_sv, col = plot_cols[5], lwd = 1.2, lty = plot_ltys[5])
lines(df9$g_k,  col = plot_cols[6], lwd = 1.2, lty = plot_ltys[6])

# --- Bottom Legend (Single Row) ---
par(mar = c(0, 0, 0, 0))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
legend("center",
       legend = plot_labels,
       col = plot_cols, lty = plot_ltys, pch = c(20, rep(NA, 5)),
       lwd = 1.2, cex =0.8, bty = "n", ncol = 6, text.width = max(strwidth(plot_labels)) * 1.0)

dev.off()