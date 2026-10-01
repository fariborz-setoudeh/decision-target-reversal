# GitHub entry point: Rscript 04_unequal_prior_simulation.R PROJECT_DIR
# Run from the repository root in RStudio, or supply its absolute path.
# ======================================================================
# Decision-Target Reversal under Informative Label Missingness
#
# SCRIPT 04: Unequal-prior robustness experiment
#
# Configuration:
#   pi1    = 0.30
#   pi2    = 0.70
#   Delta  = 1.50
#   gamma  = 0.30
#   a      = 0.25
#   xi1    = -0.25
#
# Monte Carlo:
#   n      = 2000
#   B      = 3000
#
# M = 1 means that the class label is missing.
#
# IMPORTANT:
#   - Fixed-gamma calibration is used only at the true configuration.
#   - xi0 and xi1 are freely estimated in the PC likelihood.
#   - CC and PC estimators use the same canonical Gaussian model.
#   - Exact deployment risk is evaluated analytically.
# ======================================================================


# ======================================================================
# 0. Clean workspace and packages
# ======================================================================

rm(list = ls())

library(ggplot2)
library(dplyr)
library(readr)
library(scales)


# ======================================================================
# 1. Directories
# ======================================================================

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) args[1L] else getwd()
setwd(project_dir)

out_dir <- file.path(
  getwd(),
  "decision_target_reversal_results",
  "unequal_prior_robustness"
)

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


# ======================================================================
# 2. True configuration
# ======================================================================

pi1_true <- 0.30
pi2_true <- 1 - pi1_true

Delta_true <- 1.50

lambda_true <- log(pi1_true / pi2_true)

gamma_target <- 0.30

a_true   <- 0.25
xi1_true <- -a_true

n <- 2000
B <- 3000

base_seed <- 260915


cat("\n============================================================\n")
cat("UNEQUAL-PRIOR ROBUSTNESS EXPERIMENT\n")
cat("============================================================\n")

cat("\nTrue configuration:\n")
cat("pi1     =", pi1_true, "\n")
cat("pi2     =", pi2_true, "\n")
cat("lambda  =", lambda_true, "\n")
cat("Delta   =", Delta_true, "\n")
cat("gamma   =", gamma_target, "\n")
cat("a       =", a_true, "\n")
cat("xi1     =", xi1_true, "\n")
cat("n       =", n, "\n")
cat("B       =", B, "\n")


# ======================================================================
# 3. Basic model functions
# ======================================================================

mixture_density <- function(y, pi1_value, Delta_value) {
  
  pi2_value <- 1 - pi1_value
  
  pi1_value * dnorm(
    y,
    mean = Delta_value / 2,
    sd = 1
  ) +
    pi2_value * dnorm(
      y,
      mean = -Delta_value / 2,
      sd = 1
    )
}


eta_fun <- function(y, beta0_value, Delta_value) {
  
  beta0_value + Delta_value * y
}


tau_fun <- function(y, beta0_value, Delta_value) {
  
  plogis(
    eta_fun(
      y = y,
      beta0_value = beta0_value,
      Delta_value = Delta_value
    )
  )
}


q_fun <- function(y,
                  beta0_value,
                  Delta_value,
                  xi0_value,
                  xi1_value) {
  
  eta <- eta_fun(
    y = y,
    beta0_value = beta0_value,
    Delta_value = Delta_value
  )
  
  plogis(
    xi0_value +
      xi1_value * eta^2
  )
}


# ======================================================================
# 4. Expectation under the TRUE marginal distribution of Y
#
# No recursive default arguments are used.
# ======================================================================

EY <- function(fun, rel.tol = 1e-10) {
  
  integrate(
    function(y) {
      
      fun(y) *
        mixture_density(
          y = y,
          pi1_value = pi1_true,
          Delta_value = Delta_true
        )
      
    },
    lower = -Inf,
    upper = Inf,
    subdivisions = 2000,
    rel.tol = rel.tol,
    stop.on.error = TRUE
  )$value
}


# ======================================================================
# 5. Calibrate xi0 at the TRUE configuration
#
# Solve E[q(Y)] = gamma_target.
#
# This calibration is NOT imposed during likelihood estimation.
# ======================================================================

gamma_equation <- function(xi0_value) {
  
  EY(
    function(y) {
      
      q_fun(
        y = y,
        beta0_value = lambda_true,
        Delta_value = Delta_true,
        xi0_value = xi0_value,
        xi1_value = xi1_true
      )
      
    }
  ) - gamma_target
}


xi0_true <- uniroot(
  gamma_equation,
  interval = c(-30, 30),
  tol = 1e-12
)$root


gamma_check <- EY(
  function(y) {
    
    q_fun(
      y = y,
      beta0_value = lambda_true,
      Delta_value = Delta_true,
      xi0_value = xi0_true,
      xi1_value = xi1_true
    )
    
  }
)


cat("\nCalibrated mechanism:\n")
cat("xi0     =", format(xi0_true, digits = 14), "\n")
cat("xi1     =", xi1_true, "\n")
cat("E[q(Y)] =", format(gamma_check, digits = 14), "\n")


# ======================================================================
# 6. Population information matrices
# ======================================================================

# ----------------------------------------------------------------------
# 6.1 Complete-classification information
#
# beta = (beta0, Delta)
#
# I_CC =
#
# [ pi1*pi2    0   ]
# [    0       1/4 ]
# ----------------------------------------------------------------------

I_CC <- matrix(
  c(
    pi1_true * pi2_true, 0,
    0,                    1 / 4
  ),
  nrow = 2,
  byrow = TRUE
)


# ----------------------------------------------------------------------
# 6.2 Lost conditional-label information
#
# E[q(Y) I_{Z|Y}(Y)]
# ----------------------------------------------------------------------

I_lost <- matrix(
  0,
  nrow = 2,
  ncol = 2
)


for (j in 1:2) {
  
  for (k in 1:2) {
    
    I_lost[j, k] <- EY(
      function(y) {
        
        eta <- lambda_true + Delta_true * y
        
        tau <- plogis(eta)
        
        q <- plogis(
          xi0_true +
            xi1_true * eta^2
        )
        
        xj <- if (j == 1) 1 else y
        xk <- if (k == 1) 1 else y
        
        q *
          tau *
          (1 - tau) *
          xj *
          xk
      }
    )
  }
}


# ----------------------------------------------------------------------
# 6.3 Missingness-mechanism information blocks
#
# beta-score:
#
# (M-q) 2 xi1 eta (1,y)^T
#
# nuisance xi-score:
#
# (M-q) (1, eta^2)^T
# ----------------------------------------------------------------------

I_bb_M <- matrix(
  0,
  nrow = 2,
  ncol = 2
)

I_bx_M <- matrix(
  0,
  nrow = 2,
  ncol = 2
)

I_xx_M <- matrix(
  0,
  nrow = 2,
  ncol = 2
)


for (j in 1:2) {
  
  for (k in 1:2) {
    
    # beta-beta block
    I_bb_M[j, k] <- EY(
      function(y) {
        
        eta <- lambda_true + Delta_true * y
        
        q <- plogis(
          xi0_true +
            xi1_true * eta^2
        )
        
        w <- q * (1 - q)
        
        xb_j <- if (j == 1) 1 else y
        xb_k <- if (k == 1) 1 else y
        
        gb_j <- 2 * xi1_true * eta * xb_j
        gb_k <- 2 * xi1_true * eta * xb_k
        
        w * gb_j * gb_k
      }
    )
    
    
    # beta-xi block
    I_bx_M[j, k] <- EY(
      function(y) {
        
        eta <- lambda_true + Delta_true * y
        
        q <- plogis(
          xi0_true +
            xi1_true * eta^2
        )
        
        w <- q * (1 - q)
        
        xb_j <- if (j == 1) 1 else y
        
        gb_j <- 2 * xi1_true * eta * xb_j
        
        gx_k <- if (k == 1) {
          1
        } else {
          eta^2
        }
        
        w * gb_j * gx_k
      }
    )
    
    
    # xi-xi block
    I_xx_M[j, k] <- EY(
      function(y) {
        
        eta <- lambda_true + Delta_true * y
        
        q <- plogis(
          xi0_true +
            xi1_true * eta^2
        )
        
        w <- q * (1 - q)
        
        gx_j <- if (j == 1) {
          1
        } else {
          eta^2
        }
        
        gx_k <- if (k == 1) {
          1
        } else {
          eta^2
        }
        
        w * gx_j * gx_k
      }
    )
  }
}


# ----------------------------------------------------------------------
# 6.4 Efficient missingness-mechanism information
# ----------------------------------------------------------------------

I_M_eff <- I_bb_M -
  I_bx_M %*%
  solve(I_xx_M) %*%
  t(I_bx_M)

I_M_eff <- (I_M_eff + t(I_M_eff)) / 2


# ----------------------------------------------------------------------
# 6.5 Efficient PC information
# ----------------------------------------------------------------------

I_PC <- I_CC -
  I_lost +
  I_M_eff

I_PC <- (I_PC + t(I_PC)) / 2


# ----------------------------------------------------------------------
# 6.6 Population asymptotic covariance matrices
# ----------------------------------------------------------------------

V_CC <- solve(I_CC)

V_PC <- solve(I_PC)

D_pop <- V_PC - V_CC

D_pop <- (D_pop + t(D_pop)) / 2


eig_D_pop <- sort(
  eigen(
    D_pop,
    symmetric = TRUE
  )$values
)


cat("\n============================================================\n")
cat("POPULATION INFORMATION GEOMETRY\n")
cat("============================================================\n")

cat("\nI_CC:\n")
print(I_CC, digits = 10)

cat("\nLost label information:\n")
print(I_lost, digits = 10)

cat("\nEfficient mechanism information:\n")
print(I_M_eff, digits = 10)

cat("\nI_PC:\n")
print(I_PC, digits = 10)

cat("\nV_CC:\n")
print(V_CC, digits = 10)

cat("\nV_PC:\n")
print(V_PC, digits = 10)

cat("\nD = V_PC - V_CC:\n")
print(D_pop, digits = 10)

cat("\nEigenvalues of D:\n")
print(eig_D_pop, digits = 10)


# ----------------------------------------------------------------------
# 6.7 Sanity checks
# ----------------------------------------------------------------------

if (min(eigen(I_PC, symmetric = TRUE)$values) <= 0) {
  stop("Population I_PC is not positive definite.")
}

if (!(min(eig_D_pop) < 0 && max(eig_D_pop) > 0)) {
  stop("Population D is not indefinite. Stop and diagnose.")
}


# ======================================================================
# 7. Critical-cost function
# ======================================================================

critical_costs <- function(D_matrix,
                           lambda_value,
                           Delta_value,
                           tol = 1e-10) {
  
  d00 <- D_matrix[1, 1]
  d01 <- D_matrix[1, 2]
  d11 <- D_matrix[2, 2]
  
  discriminant <- d01^2 - d00 * d11
  
  
  if (
    abs(d11) > tol &&
    discriminant > tol
  ) {
    
    t1 <- (
      -d01 -
        sqrt(discriminant)
    ) / d11
    
    t2 <- (
      -d01 +
        sqrt(discriminant)
    ) / d11
    
    roots_t <- sort(
      c(t1, t2)
    )
    
    roots_c <- exp(
      lambda_value +
        Delta_value * roots_t
    )
    
    return(
      list(
        type = "quadratic",
        discriminant = discriminant,
        t = roots_t,
        c = roots_c
      )
    )
  }
  
  
  if (
    abs(d11) <= tol &&
    abs(d01) > tol
  ) {
    
    t0 <- -d00 / (2 * d01)
    
    c0 <- exp(
      lambda_value +
        Delta_value * t0
    )
    
    return(
      list(
        type = "linear",
        discriminant = discriminant,
        t = t0,
        c = c0
      )
    )
  }
  
  
  return(
    list(
      type = "no genuine finite reversal",
      discriminant = discriminant,
      t = numeric(0),
      c = numeric(0)
    )
  )
}


pop_roots <- critical_costs(
  D_matrix = D_pop,
  lambda_value = lambda_true,
  Delta_value = Delta_true
)


cat("\nPopulation critical-cost result:\n")
print(pop_roots)


if (length(pop_roots$c) != 2) {
  stop("Expected two population critical costs but did not obtain them.")
}


c_minus_pop <- pop_roots$c[1]
c_plus_pop  <- pop_roots$c[2]


cat("\nPopulation critical costs:\n")
cat("c_minus =", format(c_minus_pop, digits = 12), "\n")
cat("c_plus  =", format(c_plus_pop, digits = 12), "\n")


# ======================================================================
# 8. Exact deployment risk
#
# C21 = 1
# C12 = c
#
# Estimated cost-sensitive boundary:
#
#   t_hat = {log(c) - beta0_hat}/Delta_hat.
#
# Risk is evaluated under the TRUE data-generating distribution.
# ======================================================================

deployment_risk <- function(beta0_hat,
                            Delta_hat,
                            c_value,
                            pi1_value,
                            Delta_true_value) {
  
  if (any(!is.finite(beta0_hat))) {
    stop("Non-finite beta0_hat in deployment_risk().")
  }
  
  if (any(!is.finite(Delta_hat))) {
    stop("Non-finite Delta_hat in deployment_risk().")
  }
  
  if (any(Delta_hat <= 0)) {
    stop("Delta_hat must be positive.")
  }
  
  if (!is.finite(c_value) || c_value <= 0) {
    stop("c_value must be positive and finite.")
  }
  
  if (
    !is.finite(pi1_value) ||
    pi1_value <= 0 ||
    pi1_value >= 1
  ) {
    stop("pi1_value must lie in (0,1).")
  }
  
  if (
    !is.finite(Delta_true_value) ||
    Delta_true_value <= 0
  ) {
    stop("Delta_true_value must be positive.")
  }
  
  
  pi2_value <- 1 - pi1_value
  
  
  threshold <- (
    log(c_value) -
      beta0_hat
  ) / Delta_hat
  
  
  # True class 1, but predict class 2
  false_negative <- pnorm(
    threshold,
    mean = Delta_true_value / 2,
    sd = 1
  )
  
  
  # True class 2, but predict class 1
  false_positive <- pnorm(
    threshold,
    mean = -Delta_true_value / 2,
    sd = 1,
    lower.tail = FALSE
  )
  
  
  # C21 = 1, C12 = c_value
  risk <- pi1_value * false_negative +
    c_value * pi2_value * false_positive
  
  
  return(risk)
}


# ======================================================================
# 9. Complete-data negative log-likelihood
#
# IMPORTANT:
# For unequal priors, beta0 and Delta are estimated jointly under the
# canonical model rather than by separate moment formulas.
#
# theta = (beta0, log Delta)
# ======================================================================

nll_CC <- function(theta, y, z) {
  
  beta0 <- theta[1]
  Delta_value <- exp(theta[2])
  
  pi1_value <- plogis(beta0)
  pi2_value <- 1 - pi1_value
  
  
  loglik <- ifelse(
    z == 1L,
    
    log(pi1_value) +
      dnorm(
        y,
        mean = Delta_value / 2,
        sd = 1,
        log = TRUE
      ),
    
    log(pi2_value) +
      dnorm(
        y,
        mean = -Delta_value / 2,
        sd = 1,
        log = TRUE
      )
  )
  
  
  value <- -sum(loglik)
  
  if (!is.finite(value)) {
    value <- 1e100
  }
  
  return(value)
}


# ======================================================================
# 10. Complete-classification MLE
# ======================================================================

fit_CC <- function(y, z) {
  
  p_start <- mean(z == 1L)
  
  p_start <- min(
    max(p_start, 1e-6),
    1 - 1e-6
  )
  
  
  s <- ifelse(
    z == 1L,
    1,
    -1
  )
  
  
  Delta_start <- 2 * mean(s * y)
  
  Delta_start <- max(
    Delta_start,
    0.10
  )
  
  
  start <- c(
    qlogis(p_start),
    log(Delta_start)
  )
  
  
  fit <- try(
    optim(
      par = start,
      fn = nll_CC,
      y = y,
      z = z,
      method = "BFGS",
      control = list(
        maxit = 500,
        reltol = 1e-11
      )
    ),
    silent = TRUE
  )
  
  
  if (
    inherits(fit, "try-error") ||
    !is.finite(fit$value)
  ) {
    
    return(
      list(
        success = FALSE
      )
    )
  }
  
  
  return(
    list(
      success = TRUE,
      beta0 = fit$par[1],
      Delta = exp(fit$par[2]),
      nll = fit$value,
      convergence = fit$convergence
    )
  )
}


# ======================================================================
# 11. Partial-classification negative log-likelihood
#
# Observed data:
#
#   O = (Y, M, (1-M)Z)
#
# Likelihood factorization:
#
# p_beta(y)
# p_beta(z|y)^(1-M)
# q_beta,xi(y)^M
# {1-q_beta,xi(y)}^(1-M)
#
# theta = (beta0, log Delta, xi0, xi1)
# ======================================================================

nll_PC <- function(theta, y, z, M) {
  
  beta0 <- theta[1]
  
  Delta_value <- exp(theta[2])
  
  xi0 <- theta[3]
  xi1 <- theta[4]
  
  
  pi1_value <- plogis(beta0)
  pi2_value <- 1 - pi1_value
  
  
  # --------------------------------------------------------------------
  # Marginal density p_beta(y)
  # --------------------------------------------------------------------
  
  logf1 <- log(pi1_value) +
    dnorm(
      y,
      mean = Delta_value / 2,
      sd = 1,
      log = TRUE
    )
  
  
  logf2 <- log(pi2_value) +
    dnorm(
      y,
      mean = -Delta_value / 2,
      sd = 1,
      log = TRUE
    )
  
  
  mx <- pmax(logf1, logf2)
  
  log_py <- mx +
    log(
      exp(logf1 - mx) +
        exp(logf2 - mx)
    )
  
  
  # --------------------------------------------------------------------
  # Posterior class probability
  # --------------------------------------------------------------------
  
  eta <- beta0 +
    Delta_value * y
  
  
  # Stable log probabilities:
  #
  # log tau       = -softplus(-eta)
  # log(1 - tau)  = -softplus( eta)
  #
  # plogis-based clipping is adequate here.
  # --------------------------------------------------------------------
  
  tau <- plogis(eta)
  
  tau <- pmin(
    pmax(tau, 1e-12),
    1 - 1e-12
  )
  
  
  # --------------------------------------------------------------------
  # Missing-label probability
  # --------------------------------------------------------------------
  
  q <- plogis(
    xi0 +
      xi1 * eta^2
  )
  
  q <- pmin(
    pmax(q, 1e-12),
    1 - 1e-12
  )
  
  
  # --------------------------------------------------------------------
  # Observed-label contribution
  # --------------------------------------------------------------------
  
  log_label <- numeric(length(y))
  
  observed <- M == 0L
  
  
  if (any(observed)) {
    
    log_label[observed] <- ifelse(
      z[observed] == 1L,
      log(tau[observed]),
      log1p(-tau[observed])
    )
  }
  
  
  # --------------------------------------------------------------------
  # Missingness-indicator contribution
  # --------------------------------------------------------------------
  
  log_M <- M * log(q) +
    (1 - M) * log1p(-q)
  
  
  ll <- log_py +
    log_label +
    log_M
  
  
  value <- -sum(ll)
  
  
  if (!is.finite(value)) {
    value <- 1e100
  }
  
  
  return(value)
}


# ======================================================================
# 12. Fit partially classified model
# ======================================================================

fit_PC <- function(y, z, M, cc_fit) {
  
  beta0_start <- cc_fit$beta0
  Delta_start <- cc_fit$Delta
  
  
  eta_start <- beta0_start +
    Delta_start * y
  
  
  # --------------------------------------------------------------------
  # Starting values for missingness regression
  # --------------------------------------------------------------------
  
  glm_fit <- try(
    suppressWarnings(
      glm(
        M ~ I(eta_start^2),
        family = binomial()
      )
    ),
    silent = TRUE
  )
  
  
  if (
    !inherits(glm_fit, "try-error") &&
    all(is.finite(coef(glm_fit)))
  ) {
    
    xi_start <- coef(glm_fit)
    
  } else {
    
    xi_start <- c(
      qlogis(
        min(
          max(mean(M), 1e-4),
          1 - 1e-4
        )
      ),
      -0.10
    )
  }
  
  
  start <- c(
    beta0_start,
    log(Delta_start),
    xi_start[1],
    xi_start[2]
  )
  
  
  # --------------------------------------------------------------------
  # Primary BFGS optimization
  # --------------------------------------------------------------------
  
  fit1 <- try(
    optim(
      par = start,
      fn = nll_PC,
      y = y,
      z = z,
      M = M,
      method = "BFGS",
      control = list(
        maxit = 800,
        reltol = 1e-10
      )
    ),
    silent = TRUE
  )
  
  
  good_fit1 <- !inherits(fit1, "try-error") &&
    is.finite(fit1$value) &&
    all(is.finite(fit1$par))
  
  
  if (good_fit1 && fit1$convergence == 0) {
    
    fit_final <- fit1
    
  } else {
    
    # ------------------------------------------------------------------
    # Fallback: Nelder-Mead
    # ------------------------------------------------------------------
    
    fit2 <- try(
      optim(
        par = start,
        fn = nll_PC,
        y = y,
        z = z,
        M = M,
        method = "Nelder-Mead",
        control = list(
          maxit = 2500,
          reltol = 1e-10
        )
      ),
      silent = TRUE
    )
    
    
    good_fit2 <- !inherits(fit2, "try-error") &&
      is.finite(fit2$value) &&
      all(is.finite(fit2$par))
    
    
    if (!good_fit2) {
      
      return(
        list(
          success = FALSE
        )
      )
    }
    
    
    # ------------------------------------------------------------------
    # BFGS polish
    # ------------------------------------------------------------------
    
    fit3 <- try(
      optim(
        par = fit2$par,
        fn = nll_PC,
        y = y,
        z = z,
        M = M,
        method = "BFGS",
        control = list(
          maxit = 800,
          reltol = 1e-10
        )
      ),
      silent = TRUE
    )
    
    
    good_fit3 <- !inherits(fit3, "try-error") &&
      is.finite(fit3$value) &&
      all(is.finite(fit3$par))
    
    
    if (
      good_fit3 &&
      fit3$value <= fit2$value + 1e-8
    ) {
      
      fit_final <- fit3
      
    } else {
      
      fit_final <- fit2
    }
  }
  
  
  pars <- fit_final$par
  
  
  return(
    list(
      success = TRUE,
      beta0 = pars[1],
      Delta = exp(pars[2]),
      xi0 = pars[3],
      xi1 = pars[4],
      nll = fit_final$value,
      convergence = fit_final$convergence
    )
  )
}


# ======================================================================
# 13. Simulate one paired replication
# ======================================================================

simulate_one <- function(b) {
  
  set.seed(
    base_seed + b
  )
  
  
  # --------------------------------------------------------------------
  # Class labels
  # --------------------------------------------------------------------
  
  z <- ifelse(
    runif(n) < pi1_true,
    1L,
    2L
  )
  
  
  # --------------------------------------------------------------------
  # Gaussian observations
  # --------------------------------------------------------------------
  
  y <- rnorm(
    n,
    mean = ifelse(
      z == 1L,
      Delta_true / 2,
      -Delta_true / 2
    ),
    sd = 1
  )
  
  
  # --------------------------------------------------------------------
  # Informative label missingness
  # --------------------------------------------------------------------
  
  eta_true <- lambda_true +
    Delta_true * y
  
  
  q_true <- plogis(
    xi0_true +
      xi1_true * eta_true^2
  )
  
  
  M <- rbinom(
    n,
    size = 1,
    prob = q_true
  )
  
  
  # --------------------------------------------------------------------
  # CC estimator
  # --------------------------------------------------------------------
  
  cc <- fit_CC(
    y = y,
    z = z
  )
  
  
  if (!cc$success) {
    
    return(
      data.frame(
        rep = b,
        success = FALSE,
        
        beta0_CC = NA_real_,
        Delta_CC = NA_real_,
        
        beta0_PC = NA_real_,
        Delta_PC = NA_real_,
        xi0_PC = NA_real_,
        xi1_PC = NA_real_,
        
        missing_fraction = mean(M),
        
        convergence_CC = NA_integer_,
        convergence_PC = NA_integer_
      )
    )
  }
  
  
  # --------------------------------------------------------------------
  # PC estimator
  # --------------------------------------------------------------------
  
  pc <- fit_PC(
    y = y,
    z = z,
    M = M,
    cc_fit = cc
  )
  
  
  if (!pc$success) {
    
    return(
      data.frame(
        rep = b,
        success = FALSE,
        
        beta0_CC = cc$beta0,
        Delta_CC = cc$Delta,
        
        beta0_PC = NA_real_,
        Delta_PC = NA_real_,
        xi0_PC = NA_real_,
        xi1_PC = NA_real_,
        
        missing_fraction = mean(M),
        
        convergence_CC = cc$convergence,
        convergence_PC = NA_integer_
      )
    )
  }
  
  
  return(
    data.frame(
      rep = b,
      success = TRUE,
      
      beta0_CC = cc$beta0,
      Delta_CC = cc$Delta,
      
      beta0_PC = pc$beta0,
      Delta_PC = pc$Delta,
      
      xi0_PC = pc$xi0,
      xi1_PC = pc$xi1,
      
      missing_fraction = mean(M),
      
      convergence_CC = cc$convergence,
      convergence_PC = pc$convergence
    )
  )
}


# ======================================================================
# 14. Monte Carlo simulation with checkpointing
#
# IMPORTANT:
# This corrected version uses a new checkpoint filename so that an
# earlier run based on the old CC estimator cannot be accidentally reused.
# ======================================================================

checkpoint_file <- file.path(
  out_dir,
  "unequal_prior_checkpoint_CORRECTED.rds"
)


if (file.exists(checkpoint_file)) {
  
  sim_results <- readRDS(
    checkpoint_file
  )
  
  
  if (nrow(sim_results) > 0) {
    
    start_rep <- max(sim_results$rep) + 1
    
  } else {
    
    start_rep <- 1
  }
  
  
  cat(
    "\nCorrected checkpoint found. Resuming at replication",
    start_rep,
    "\n"
  )
  
} else {
  
  sim_results <- data.frame()
  
  start_rep <- 1
}


if (start_rep <= B) {
  
  for (b in start_rep:B) {
    
    one_result <- simulate_one(b)
    
    
    sim_results <- bind_rows(
      sim_results,
      one_result
    )
    
    
    if (
      b %% 25 == 0 ||
      b == B
    ) {
      
      cat(
        sprintf(
          "\rReplication %4d / %4d | successful = %4d",
          b,
          B,
          sum(sim_results$success)
        )
      )
      
      flush.console()
    }
    
    
    if (
      b %% 100 == 0 ||
      b == B
    ) {
      
      saveRDS(
        sim_results,
        checkpoint_file
      )
    }
  }
}


cat("\n")


# ======================================================================
# 15. Save raw corrected results
# ======================================================================

sim_results <- sim_results %>%
  arrange(rep)


write_csv(
  sim_results,
  file.path(
    out_dir,
    "unequal_prior_raw_results_CORRECTED.csv"
  )
)


saveRDS(
  sim_results,
  file.path(
    out_dir,
    "unequal_prior_raw_results_CORRECTED.rds"
  )
)


# ======================================================================
# 16. Successful replications
# ======================================================================

valid <- sim_results %>%
  filter(success)


B_success <- nrow(valid)

success_rate <- B_success / B


cat("\n============================================================\n")
cat("MONTE CARLO DIAGNOSTICS\n")
cat("============================================================\n")

cat("\nSuccessful replications =", B_success, "/", B, "\n")
cat("Success rate            =", success_rate, "\n")

cat(
  "Mean missing fraction  =",
  format(
    mean(valid$missing_fraction),
    digits = 10
  ),
  "\n"
)


if (B_success < 0.95 * B) {
  warning("Less than 95% of the replications were successfully fitted.")
}


# ======================================================================
# 17. Parameter summaries
# ======================================================================

parameter_summary <- data.frame(
  
  parameter = c(
    "beta0_CC",
    "Delta_CC",
    "beta0_PC",
    "Delta_PC",
    "xi0_PC",
    "xi1_PC"
  ),
  
  truth = c(
    lambda_true,
    Delta_true,
    lambda_true,
    Delta_true,
    xi0_true,
    xi1_true
  ),
  
  mean = c(
    mean(valid$beta0_CC),
    mean(valid$Delta_CC),
    mean(valid$beta0_PC),
    mean(valid$Delta_PC),
    mean(valid$xi0_PC),
    mean(valid$xi1_PC)
  ),
  
  sd = c(
    sd(valid$beta0_CC),
    sd(valid$Delta_CC),
    sd(valid$beta0_PC),
    sd(valid$Delta_PC),
    sd(valid$xi0_PC),
    sd(valid$xi1_PC)
  )
)


parameter_summary$bias <- parameter_summary$mean -
  parameter_summary$truth


cat("\nParameter summary:\n")

print(
  parameter_summary,
  digits = 8,
  row.names = FALSE
)


write_csv(
  parameter_summary,
  file.path(
    out_dir,
    "unequal_prior_parameter_summary_CORRECTED.csv"
  )
)


# ======================================================================
# 18. Empirical covariance geometry
#
# V is covariance of sqrt(n)(beta_hat-beta0), hence n*Cov(beta_hat).
# ======================================================================

beta_CC <- cbind(
  valid$beta0_CC,
  valid$Delta_CC
)


beta_PC <- cbind(
  valid$beta0_PC,
  valid$Delta_PC
)


V_CC_emp <- n * cov(beta_CC)

V_PC_emp <- n * cov(beta_PC)


D_emp <- V_PC_emp -
  V_CC_emp

D_emp <- (D_emp + t(D_emp)) / 2


eig_D_emp <- sort(
  eigen(
    D_emp,
    symmetric = TRUE
  )$values
)


cat("\n============================================================\n")
cat("EMPIRICAL COVARIANCE GEOMETRY\n")
cat("============================================================\n")

cat("\nEmpirical V_CC:\n")
print(V_CC_emp, digits = 10)

cat("\nEmpirical V_PC:\n")
print(V_PC_emp, digits = 10)

cat("\nEmpirical D:\n")
print(D_emp, digits = 10)

cat("\nEigenvalues of empirical D:\n")
print(eig_D_emp, digits = 10)


emp_roots <- critical_costs(
  D_matrix = D_emp,
  lambda_value = lambda_true,
  Delta_value = Delta_true
)


cat("\nEmpirical-D critical-cost result:\n")
print(emp_roots)


# ======================================================================
# 19. Cost grid
# ======================================================================

log_c_grid <- seq(
  from = log(c_minus_pop) - 1,
  to = log(c_plus_pop) + 1,
  length.out = 181
)


c_grid <- exp(log_c_grid)


# ======================================================================
# 20. Exact finite-sample risk curves and first-order approximations
# ======================================================================

risk_summary_list <- vector(
  "list",
  length(c_grid)
)


for (g in seq_along(c_grid)) {
  
  c_now <- c_grid[g]
  
  
  # --------------------------------------------------------------------
  # Exact deployment risks
  # --------------------------------------------------------------------
  
  risk_CC <- deployment_risk(
    beta0_hat = valid$beta0_CC,
    Delta_hat = valid$Delta_CC,
    c_value = c_now,
    pi1_value = pi1_true,
    Delta_true_value = Delta_true
  )
  
  
  risk_PC <- deployment_risk(
    beta0_hat = valid$beta0_PC,
    Delta_hat = valid$Delta_PC,
    c_value = c_now,
    pi1_value = pi1_true,
    Delta_true_value = Delta_true
  )
  
  
  risk_diff <- risk_PC -
    risk_CC
  
  
  mean_diff <- mean(risk_diff)
  
  se_diff <- sd(risk_diff) /
    sqrt(length(risk_diff))
  
  
  # --------------------------------------------------------------------
  # True cost-sensitive Bayes boundary
  # --------------------------------------------------------------------
  
  tC <- (
    log(c_now) -
      lambda_true
  ) / Delta_true
  
  
  # --------------------------------------------------------------------
  # Boundary weight
  #
  # C21 = 1:
  #
  # m_C = pi1 phi(t_C - Delta/2)
  #
  # This equals c*pi2*phi(t_C+Delta/2) at the Bayes boundary.
  # --------------------------------------------------------------------
  
  mC_left <- pi1_true *
    dnorm(
      tC -
        Delta_true / 2
    )
  
  
  mC_right <- c_now *
    pi2_true *
    dnorm(
      tC +
        Delta_true / 2
    )
  
  
  # Numerical identity check
  if (
    abs(mC_left - mC_right) >
    1e-8 * max(1, abs(mC_left), abs(mC_right))
  ) {
    
    stop(
      "Boundary identity for m_C failed."
    )
  }
  
  
  mC <- mC_left
  
  
  # --------------------------------------------------------------------
  # Population quadratic form
  # --------------------------------------------------------------------
  
  Q_pop <- D_pop[1, 1] +
    2 * tC * D_pop[1, 2] +
    tC^2 * D_pop[2, 2]
  
  
  population_theory <- (
    mC /
      (2 * Delta_true * n)
  ) * Q_pop
  
  
  # --------------------------------------------------------------------
  # Empirical-D quadratic form
  # --------------------------------------------------------------------
  
  Q_emp <- D_emp[1, 1] +
    2 * tC * D_emp[1, 2] +
    tC^2 * D_emp[2, 2]
  
  
  empirical_D_theory <- (
    mC /
      (2 * Delta_true * n)
  ) * Q_emp
  
  
  risk_summary_list[[g]] <- data.frame(
    
    c = c_now,
    log_c = log(c_now),
    
    mean_diff = mean_diff,
    se_diff = se_diff,
    
    lower_95 = mean_diff -
      1.96 * se_diff,
    
    upper_95 = mean_diff +
      1.96 * se_diff,
    
    population_theory = population_theory,
    
    empirical_D_theory = empirical_D_theory,
    
    n_mean_diff = n * mean_diff,
    
    n_population_theory = n * population_theory,
    
    n_empirical_D_theory = n * empirical_D_theory
  )
}


risk_summary <- bind_rows(
  risk_summary_list
)


write_csv(
  risk_summary,
  file.path(
    out_dir,
    "unequal_prior_risk_summary_CORRECTED.csv"
  )
)


# ======================================================================
# 21. Locate zero crossings of Monte Carlo risk curve
# ======================================================================

find_crossings <- function(x, y) {
  
  roots <- numeric(0)
  
  
  for (i in seq_len(length(y) - 1)) {
    
    if (
      is.finite(y[i]) &&
      is.finite(y[i + 1])
    ) {
      
      if (y[i] == 0) {
        
        roots <- c(
          roots,
          x[i]
        )
        
      } else if (
        y[i] * y[i + 1] < 0
      ) {
        
        x1 <- x[i]
        x2 <- x[i + 1]
        
        y1 <- y[i]
        y2 <- y[i + 1]
        
        
        root <- x1 -
          y1 *
          (x2 - x1) /
          (y2 - y1)
        
        
        roots <- c(
          roots,
          root
        )
      }
    }
  }
  
  
  return(
    exp(
      unique(roots)
    )
  )
}


MC_roots <- find_crossings(
  x = risk_summary$log_c,
  y = risk_summary$mean_diff
)


cat("\n============================================================\n")
cat("CRITICAL-COST COMPARISON\n")
cat("============================================================\n")

cat("\nPopulation roots:\n")
print(pop_roots$c, digits = 10)

cat("\nEmpirical-D roots:\n")
print(emp_roots$c, digits = 10)

cat("\nMonte Carlo risk-curve crossings:\n")
print(MC_roots, digits = 10)


# ======================================================================
# 22. Equal-cost comparison
# ======================================================================

risk_CC_equal <- deployment_risk(
  beta0_hat = valid$beta0_CC,
  Delta_hat = valid$Delta_CC,
  c_value = 1,
  pi1_value = pi1_true,
  Delta_true_value = Delta_true
)


risk_PC_equal <- deployment_risk(
  beta0_hat = valid$beta0_PC,
  Delta_hat = valid$Delta_PC,
  c_value = 1,
  pi1_value = pi1_true,
  Delta_true_value = Delta_true
)


diff_equal <- risk_PC_equal -
  risk_CC_equal


mean_equal <- mean(diff_equal)


se_equal <- sd(diff_equal) /
  sqrt(length(diff_equal))


CI_equal <- c(
  mean_equal -
    1.96 * se_equal,
  
  mean_equal +
    1.96 * se_equal
)


# Population first-order value at c = 1

t_equal <- (
  0 -
    lambda_true
) / Delta_true


m_equal <- pi1_true *
  dnorm(
    t_equal -
      Delta_true / 2
  )


Q_equal_pop <- D_pop[1, 1] +
  2 * t_equal * D_pop[1, 2] +
  t_equal^2 * D_pop[2, 2]


equal_population_theory <- (
  m_equal /
    (2 * Delta_true * n)
) * Q_equal_pop


Q_equal_emp <- D_emp[1, 1] +
  2 * t_equal * D_emp[1, 2] +
  t_equal^2 * D_emp[2, 2]


equal_empirical_theory <- (
  m_equal /
    (2 * Delta_true * n)
) * Q_equal_emp


cat("\nEqual-cost result, c = 1:\n")

cat(
  "Monte Carlo mean PC - CC =",
  format(mean_equal, scientific = TRUE, digits = 8),
  "\n"
)

cat(
  "95% MC interval = (",
  format(CI_equal[1], scientific = TRUE, digits = 8),
  ", ",
  format(CI_equal[2], scientific = TRUE, digits = 8),
  ")\n",
  sep = ""
)

cat(
  "Population first-order =",
  format(equal_population_theory, scientific = TRUE, digits = 8),
  "\n"
)

cat(
  "Empirical-D first-order =",
  format(equal_empirical_theory, scientific = TRUE, digits = 8),
  "\n"
)


# ======================================================================
# 23. RMSE of first-order approximations
# ======================================================================

RMSE_population <- sqrt(
  mean(
    (
      risk_summary$mean_diff -
        risk_summary$population_theory
    )^2
  )
)


RMSE_empiricalD <- sqrt(
  mean(
    (
      risk_summary$mean_diff -
        risk_summary$empirical_D_theory
    )^2
  )
)


cat("\nRisk-curve approximation RMSE:\n")

cat(
  "Population D =",
  format(
    RMSE_population,
    scientific = TRUE,
    digits = 8
  ),
  "\n"
)

cat(
  "Empirical D  =",
  format(
    RMSE_empiricalD,
    scientific = TRUE,
    digits = 8
  ),
  "\n"
)


# ======================================================================
# 24. Compact final summary
# ======================================================================

get_root <- function(root_object, index) {
  
  if (
    length(root_object$c) >= index &&
    is.finite(root_object$c[index])
  ) {
    
    return(root_object$c[index])
  }
  
  return(NA_real_)
}


get_MC_root <- function(root_vector, index) {
  
  if (
    length(root_vector) >= index &&
    is.finite(root_vector[index])
  ) {
    
    return(root_vector[index])
  }
  
  return(NA_real_)
}


summary_results <- data.frame(
  
  quantity = c(
    "pi1",
    "pi2",
    "lambda",
    "Delta",
    "gamma",
    "a",
    "xi0",
    "xi1",
    
    "n",
    "B",
    "successful_B",
    "success_rate",
    "mean_missing_fraction",
    
    "population_D_00",
    "population_D_01",
    "population_D_11",
    
    "population_eigenvalue_min",
    "population_eigenvalue_max",
    
    "empirical_D_00",
    "empirical_D_01",
    "empirical_D_11",
    
    "empirical_eigenvalue_min",
    "empirical_eigenvalue_max",
    
    "population_c_minus",
    "population_c_plus",
    
    "empirical_c_minus",
    "empirical_c_plus",
    
    "MC_c_minus",
    "MC_c_plus",
    
    "equal_cost_MC_difference",
    "equal_cost_CI_lower",
    "equal_cost_CI_upper",
    
    "equal_cost_population_theory",
    "equal_cost_empirical_D_theory",
    
    "RMSE_population_D",
    "RMSE_empirical_D"
  ),
  
  
  value = c(
    pi1_true,
    pi2_true,
    lambda_true,
    Delta_true,
    gamma_target,
    a_true,
    xi0_true,
    xi1_true,
    
    n,
    B,
    B_success,
    success_rate,
    mean(valid$missing_fraction),
    
    D_pop[1, 1],
    D_pop[1, 2],
    D_pop[2, 2],
    
    min(eig_D_pop),
    max(eig_D_pop),
    
    D_emp[1, 1],
    D_emp[1, 2],
    D_emp[2, 2],
    
    min(eig_D_emp),
    max(eig_D_emp),
    
    get_root(pop_roots, 1),
    get_root(pop_roots, 2),
    
    get_root(emp_roots, 1),
    get_root(emp_roots, 2),
    
    get_MC_root(MC_roots, 1),
    get_MC_root(MC_roots, 2),
    
    mean_equal,
    CI_equal[1],
    CI_equal[2],
    
    equal_population_theory,
    equal_empirical_theory,
    
    RMSE_population,
    RMSE_empiricalD
  )
)


write_csv(
  summary_results,
  file.path(
    out_dir,
    "unequal_prior_final_summary_CORRECTED.csv"
  )
)


# ======================================================================
# 25. Professional supplementary figure
# ======================================================================

plot_long <- bind_rows(
  
  data.frame(
    log_c = risk_summary$log_c,
    value = risk_summary$n_mean_diff,
    curve = "Monte Carlo"
  ),
  
  data.frame(
    log_c = risk_summary$log_c,
    value = risk_summary$n_empirical_D_theory,
    curve = "First-order using empirical covariance"
  ),
  
  data.frame(
    log_c = risk_summary$log_c,
    value = risk_summary$n_population_theory,
    curve = "Population first-order theory"
  )
)


plot_long$curve <- factor(
  plot_long$curve,
  levels = c(
    "Monte Carlo",
    "First-order using empirical covariance",
    "Population first-order theory"
  )
)


curve_colors <- c(
  "Monte Carlo" = "#0072B2",
  "First-order using empirical covariance" = "#D55E00",
  "Population first-order theory" = "#009E73"
)


curve_types <- c(
  "Monte Carlo" = "solid",
  "First-order using empirical covariance" = "longdash",
  "Population first-order theory" = "dotted"
)


curve_widths <- c(
  "Monte Carlo" = 1.25,
  "First-order using empirical covariance" = 1.15,
  "Population first-order theory" = 1.15
)


y_min <- min(
  plot_long$value,
  na.rm = TRUE
)

y_max <- max(
  plot_long$value,
  na.rm = TRUE
)

y_range <- y_max - y_min

if (!is.finite(y_range) || y_range <= 0) {
  y_range <- 1
}


shade_mid <- mean(
  c(
    log(c_minus_pop),
    log(c_plus_pop)
  )
)


p <- ggplot(
  plot_long,
  aes(
    x = log_c,
    y = value,
    color = curve,
    linetype = curve,
    linewidth = curve
  )
) +
  
  # Population-theory favourable interval
  annotate(
    "rect",
    xmin = log(c_minus_pop),
    xmax = log(c_plus_pop),
    ymin = -Inf,
    ymax = Inf,
    fill = "#56B4E9",
    alpha = 0.08
  ) +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.50,
    linetype = "dashed",
    color = "grey35"
  ) +
  
  geom_vline(
    xintercept = log(c_minus_pop),
    linewidth = 0.55,
    linetype = "dashed",
    color = "grey40"
  ) +
  
  geom_vline(
    xintercept = log(c_plus_pop),
    linewidth = 0.55,
    linetype = "dashed",
    color = "grey40"
  ) +
  
  geom_line(
    lineend = "round"
  ) +
  
  scale_color_manual(
    values = curve_colors,
    name = NULL
  ) +
  
  scale_linetype_manual(
    values = curve_types,
    name = NULL
  ) +
  
  scale_linewidth_manual(
    values = curve_widths,
    name = NULL
  ) +
  
  labs(
    x = expression(
      log(C[12] / C[21])
    ),
    
    y = expression(
      n *
        {
          E[
            R[C](hat(beta)[PC])
          ] -
            E[
              R[C](hat(beta)[CC])
            ]
        }
    )
  ) +
  
  annotate(
    "text",
    x = shade_mid,
    y = y_min + 0.10 * y_range,
    label = "Population-theory favourable region",
    color = "#0072B2",
    size = 4.2,
    fontface = "italic"
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    axis.title.x = element_text(
      size = 14,
      margin = margin(t = 10)
    ),
    
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    
    axis.text = element_text(
      size = 11.5,
      color = "black"
    ),
    
    axis.line = element_line(
      linewidth = 0.50,
      color = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.45,
      color = "black"
    ),
    
    legend.position = c(0.72, 0.86),
    
    legend.background = element_rect(
      fill = alpha("white", 0.90),
      color = NA
    ),
    
    legend.text = element_text(
      size = 11
    ),
    
    plot.margin = margin(
      t = 12,
      r = 18,
      b = 10,
      l = 10
    )
  ) +
  
  guides(
    linewidth = "none",
    
    color = guide_legend(
      override.aes = list(
        linewidth = 1.4
      )
    )
  )


print(p)


# ======================================================================
# 26. Save figure
# ======================================================================

png_file <- file.path(
  out_dir,
  "Figure_Unequal_Prior_Robustness_CORRECTED.png"
)


pdf_file <- file.path(
  out_dir,
  "Figure_Unequal_Prior_Robustness_CORRECTED.pdf"
)


ggsave(
  filename = png_file,
  plot = p,
  width = 8.4,
  height = 6.2,
  dpi = 600,
  bg = "white"
)


ggsave(
  filename = pdf_file,
  plot = p,
  width = 8.4,
  height = 6.2,
  device = cairo_pdf,
  bg = "white"
)


# ======================================================================
# 27. Save session information
# ======================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    out_dir,
    "sessionInfo_CORRECTED.txt"
  )
)


# ======================================================================
# 28. Final report
# ======================================================================

cat("\n\n============================================================\n")
cat("UNEQUAL-PRIOR ROBUSTNESS STUDY COMPLETE\n")
cat("============================================================\n")


cat("\nTRUE MECHANISM\n")
cat("xi0 =", format(xi0_true, digits = 12), "\n")
cat("xi1 =", xi1_true, "\n")
cat("E[q(Y)] =", format(gamma_check, digits = 12), "\n")


cat("\nPOPULATION D\n")
print(D_pop, digits = 10)


cat("\nPOPULATION EIGENVALUES\n")
print(eig_D_pop, digits = 10)


cat("\nPOPULATION CRITICAL COSTS\n")
print(pop_roots$c, digits = 10)


cat("\nMONTE CARLO DIAGNOSTICS\n")
cat("Successful =", B_success, "/", B, "\n")
cat(
  "Missing fraction =",
  format(mean(valid$missing_fraction), digits = 10),
  "\n"
)


cat("\nEMPIRICAL D\n")
print(D_emp, digits = 10)


cat("\nEMPIRICAL EIGENVALUES\n")
print(eig_D_emp, digits = 10)


cat("\nEMPIRICAL-D CRITICAL COSTS\n")
print(emp_roots$c, digits = 10)


cat("\nMONTE CARLO CROSSINGS\n")
print(MC_roots, digits = 10)


cat("\nEQUAL-COST RESULT\n")

cat(
  "MC difference =",
  format(mean_equal, scientific = TRUE, digits = 8),
  "\n"
)

cat(
  "95% CI = (",
  format(CI_equal[1], scientific = TRUE, digits = 8),
  ", ",
  format(CI_equal[2], scientific = TRUE, digits = 8),
  ")\n",
  sep = ""
)

cat(
  "Population theory =",
  format(equal_population_theory, scientific = TRUE, digits = 8),
  "\n"
)

cat(
  "Empirical-D theory =",
  format(equal_empirical_theory, scientific = TRUE, digits = 8),
  "\n"
)


cat("\nRISK-CURVE RMSE\n")

cat(
  "Population D =",
  format(RMSE_population, scientific = TRUE, digits = 8),
  "\n"
)

cat(
  "Empirical D  =",
  format(RMSE_empiricalD, scientific = TRUE, digits = 8),
  "\n"
)


cat("\nOUTPUT DIRECTORY\n")
cat(out_dir, "\n")


cat("\nFIGURE\n")
cat(pdf_file, "\n")


cat("\n============================================================\n")
cat("END\n")
cat("============================================================\n")