# FSSM: State-Space Modeling of Shape-constrained Functional Time Series

This repository provides the R implementation for the simulation studies proposed in the following paper. 

Hiraki, D., Hamura, Y., Irie, K. and Sugasawa, S. (2024). State-Space Modeling of Shape-constrained Functional Time Series. [arXiv:2404.07586](https://arxiv.org/abs/2404.07586)

This repository contains R scripts to implement the functional state-space model (FSSM) under monotonicity and convexity constraints. It is designed to replicate the one-shot simulation experiments ($K=4, 9$) and evaluate Gini coefficients and parameter estimations against various benchmarking models.

## File Structure

### Core Model Implementations (MCMC Samplers)
- `FSSM.R`: Implements the proposed shape-constrained functional state-space model via data augmentation and scalar FFBS.
- `Mixture.R`: Implements the mixture extension of FSSM (Basis 1) with time-varying/component-specific variances.
- `ARGP-function.R`: Implements the Functional Time Series with Gaussian Process AR(1) (ARGP) model.
- `RWIW.R`: Implements the Random Walk Inverse Wishart (RWIW) Dynamic Linear Model.
- `RWSV.R`: Implements the Random Walk Stochastic Volatility (RWSV) Dynamic Linear Model.
- `tsdir_sm.R`: Implements Kobayashi's Time-Series Dirichlet (TS-DIR) model for smooth Lorenz curves.

### Execution & Replication Script
- `main_sim_oneshot.R`: Executes the single-run simulation study to compare Gini coefficients and parameter estimations across all models for $K=4$ and $K=9$.
