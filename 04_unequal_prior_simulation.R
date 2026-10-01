
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

required_packages <- c("ggplot2", "dplyr", "readr", "scales")
missing_packages <- required_packages[!vapply(required_packages,
  requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install these packages first: ", paste(missing_packages, collapse = ", "))
}
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

run_version <- "observed-label-only-v3-20261001"
manifest_file <- file.path(out_dir, "observed_label_only_manifest.rds")
if (dir.exists(out_dir) && !file.exists(manifest_file)) {
  archive <- paste0(out_dir, "_BEFORE_OBSERVED_LABEL_FIX_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  while (file.exists(archive)) archive <- paste0(archive, "_previous")
  if (!file.rename(out_dir, archive)) {
    stop("Could not archive old results. Close files in that folder and rerun.")
  }
  cat("Previous results retained at:\n", archive, "\n", sep = "")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (file.exists(manifest_file)) {
  if (!identical(readRDS(manifest_file), run_version)) {
    stop("Output folder belongs to another script version. Use a fresh project folder.")
  }
} else {
  saveRDS(run_version, manifest_file)
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
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

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
# 9. Exact CC MLE and observed-data-only PC fitting
# Both b and Delta are estimated; midpoint and within-class variance are fixed.
# Missingness parameters remain free, without a sign restriction or penalty.

fit_CC <- function(y, z) {
  pi_hat <- mean(z == 1L)
  Delta_hat <- 2 * mean(ifelse(z == 1L, 1, -1) * y)
  if (!is.finite(pi_hat) || pi_hat <= 0 || pi_hat >= 1 ||
      !is.finite(Delta_hat) || Delta_hat <= 0) {
    return(list(success = FALSE))
  }
  b_hat <- qlogis(pi_hat)
  value <- -sum(ifelse(z == 1L, log(pi_hat), log1p(-pi_hat)) +
    dnorm(y, mean = ifelse(z == 1L, Delta_hat / 2, -Delta_hat / 2),
          sd = 1, log = TRUE))
  list(success = TRUE, beta0 = b_hat, Delta = Delta_hat,
       nll = value, convergence = 0L)
}

softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))

# Returns the negative MEAN observed log likelihood and its analytic gradient.
# Values of z at M=1 are never used; the caller replaces them by NA explicitly.
pc_value_gradient <- function(theta, y, z, M) {
  bad <- list(value = 1e100, gradient = rep(0, 4))
  if (length(theta) != 4L || !all(is.finite(theta))) return(bad)
  b <- theta[1]; Delta <- exp(theta[2])
  x0 <- theta[3]; x1 <- theta[4]
  if (!is.finite(Delta) || Delta <= 0) return(bad)
  eta <- b + Delta * y
  nu <- x0 + x1 * eta^2
  if (!all(is.finite(eta)) || !all(is.finite(nu))) return(bad)
  lp1 <- -softplus(-b) + dnorm(y, Delta / 2, 1, log = TRUE)
  lp2 <- -softplus(b) + dnorm(y, -Delta / 2, 1, log = TRUE)
  mx <- pmax(lp1, lp2)
  log_y <- mx + log(exp(lp1 - mx) + exp(lp2 - mx))
  tau <- plogis(eta); q <- plogis(nu); pi_hat <- plogis(b)
  observed <- M == 0L
  log_label <- numeric(length(y)); label_score <- numeric(length(y))
  if (any(observed)) {
    z1 <- as.numeric(z[observed] == 1L)
    log_label[observed] <- ifelse(z1 == 1, -softplus(-eta[observed]),
                                -softplus(eta[observed]))
    label_score[observed] <- z1 - tau[observed]
  }
  log_M <- ifelse(M == 1L, -softplus(-nu), -softplus(nu))
  residual_M <- M - q
  mechanism_score <- 2 * x1 * eta * residual_M
  sb <- tau - pi_hat + label_score + mechanism_score
  sd <- y * (2 * tau - 1) / 2 - Delta / 4 +
        y * (label_score + mechanism_score)
  gradient <- -c(mean(sb), Delta * mean(sd), mean(residual_M),
                 mean(residual_M * eta^2))
  value <- -mean(log_y + log_label + log_M)
  if (!is.finite(value) || !all(is.finite(gradient))) return(bad)
  list(value = value, gradient = gradient)
}

nll_PC <- function(theta, y, z, M) pc_value_gradient(theta, y, z, M)$value
grad_PC <- function(theta, y, z, M) pc_value_gradient(theta, y, z, M)$gradient

pc_starts <- function(y, z, M) {
  observed <- M == 0L
  zo <- z[observed]
  pi_start <- if (length(zo)) mean(zo == 1L) else 0.5
  pi_start <- min(max(pi_start, 1e-4), 1 - 1e-4)
  b_start <- qlogis(pi_start)
  Delta_start <- if (length(zo)) {
    2 * mean(ifelse(zo == 1L, 1, -1) * y[observed])
  } else 1
  if (!is.finite(Delta_start) || Delta_start <= 0.20) Delta_start <- 1
  Delta_start <- min(Delta_start, 5)
  eta_start <- b_start + Delta_start * y
  intercept <- qlogis(min(max(mean(M), 1e-4), 1 - 1e-4))
  glm_fit <- try(suppressWarnings(glm(M ~ I(eta_start^2), family = binomial())),
                 silent = TRUE)
  xi <- c(intercept, -0.10)
  if (!inherits(glm_fit, 'try-error') && length(coef(glm_fit)) == 2L &&
      all(is.finite(coef(glm_fit)))) xi <- unname(coef(glm_fit))
  # These bounds apply to starting values only, not to fitted parameters.
  xi <- c(min(max(xi[1], -10), 10), min(max(xi[2], -3), 1))
  rbind(c(b_start, log(Delta_start), xi),
        c(b_start, log(Delta_start), intercept, -0.10),
        c(b_start, log(Delta_start), intercept, -0.30))
}

# Newton polishing improves stationarity after a successful optimizer return.
# Every step must reduce the unpenalized observed-data objective.
polish_PC <- function(theta, y, z, M) {
  for (iteration in seq_len(6L)) {
    vg <- pc_value_gradient(theta, y, z, M)
    if (max(abs(vg$gradient)) <= 1e-9) break
    H <- try(optimHess(theta, nll_PC, gr = grad_PC, y = y, z = z, M = M),
             silent = TRUE)
    if (inherits(H, 'try-error') || !all(is.finite(H))) break
    H <- (H + t(H)) / 2
    if (min(eigen(H, symmetric = TRUE, only.values = TRUE)$values) <= 1e-8) break
    step <- try(solve(H, vg$gradient), silent = TRUE)
    if (inherits(step, 'try-error') || !all(is.finite(step))) break
    accepted <- FALSE
    for (k in 0:15) {
      alpha <- 2^(-k)
      proposal <- theta - alpha * step
      value <- nll_PC(proposal, y, z, M)
      if (is.finite(value) &&
          value <= vg$value - 1e-4 * alpha * sum(vg$gradient * step)) {
        theta <- proposal; accepted <- TRUE; break
      }
    }
    if (!accepted) break
  }
  theta
}

fit_PC <- function(y, z, M) {
  stopifnot(length(y) == length(z), length(y) == length(M),
            all(is.finite(y)), all(M %in% c(0L, 1L)),
            all(z[M == 0L] %in% c(1L, 2L)), all(is.na(z[M == 1L])))
  starts <- pc_starts(y, z, M)
  candidates <- vector('list', nrow(starts))
  for (j in seq_len(nrow(starts))) {
    start <- starts[j, ]
    opt <- try(optim(start, nll_PC, gr = grad_PC, y = y, z = z, M = M,
                     method = 'BFGS', control = list(maxit = 1500, reltol = 1e-12)),
               silent = TRUE)
    usable <- !inherits(opt, 'try-error') && is.finite(opt$value) &&
              all(is.finite(opt$par)) && opt$convergence == 0L
    if (!usable) {
      alternate <- try(nlminb(start, nll_PC, gradient = grad_PC, y = y, z = z,
                              M = M, control = list(iter.max = 1500,
                              eval.max = 4000, rel.tol = 1e-12, x.tol = 1e-10)),
                        silent = TRUE)
      if (!inherits(alternate, 'try-error') && is.finite(alternate$objective) &&
          all(is.finite(alternate$par)) && alternate$convergence == 0L) {
        opt <- list(par = alternate$par, value = alternate$objective,
                    convergence = alternate$convergence)
        usable <- TRUE
      }
    }
    if (!usable) {
      candidates[[j]] <- list(eligible = FALSE, value = Inf)
      next
    }
    theta <- polish_PC(opt$par, y, z, M)
    vg <- pc_value_gradient(theta, y, z, M)
    H <- try(optimHess(theta, nll_PC, gr = grad_PC, y = y, z = z, M = M,
                       control = list(ndeps = rep(1e-4, 4))), silent = TRUE)
    min_eigen <- if (!inherits(H, 'try-error') && all(is.finite(H))) {
      min(eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values)
    } else NA_real_
    max_gradient <- max(abs(vg$gradient))
    eligible <- is.finite(vg$value) && vg$value < 1e99 &&
                is.finite(min_eigen) && min_eigen > 1e-8 && max_gradient <= 1e-7
    candidates[[j]] <- list(eligible = eligible, value = vg$value, theta = theta,
                            convergence = opt$convergence,
                            max_gradient = max_gradient, min_eigen = min_eigen)
  }
  eligible <- vapply(candidates, function(x) isTRUE(x$eligible), logical(1))
  if (!any(eligible)) return(list(success = FALSE))
  values <- vapply(candidates, function(x) x$value, numeric(1))
  ids <- which(eligible); best_id <- ids[which.min(values[ids])]
  best <- candidates[[best_id]]
  reproduced <- sum(eligible & abs(values - best$value) <= 1e-8)
  unresolved_better <- any(!eligible & is.finite(values) & values < best$value - 1e-8)
  if (reproduced < 2L || unresolved_better) return(list(success = FALSE))
  theta <- best$theta
  list(success = TRUE, beta0 = theta[1], Delta = exp(theta[2]),
       xi0 = theta[3], xi1 = theta[4], nll = best$value * length(y),
       convergence = best$convergence, max_gradient = best$max_gradient,
       min_eigen = best$min_eigen, reproduced_starts = reproduced)
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
        convergence_PC = NA_integer_,
        pc_max_abs_gradient = NA_real_, pc_min_hessian_eigenvalue = NA_real_,
        pc_reproduced_starts = NA_integer_
      )
    )
  }

  # --------------------------------------------------------------------
  # PC estimator
  # --------------------------------------------------------------------
  
  z_PC <- z
  z_PC[M == 1L] <- NA_integer_

  if (b == 1L) {
    initial <- pc_starts(y, z_PC, M)[1, ]
    numerical <- vapply(seq_len(4L), function(j) {
      shift <- rep(0, 4); shift[j] <- 1e-5
      (nll_PC(initial + shift, y, z_PC, M) -
       nll_PC(initial - shift, y, z_PC, M)) / 2e-5
    }, numeric(1))
    analytic <- grad_PC(initial, y, z_PC, M)
    gradient_error <- max(abs(numerical - analytic)) / (1 + max(abs(analytic)))
    z_altered <- z_PC; z_altered[M == 1L] <- 2L
    stopifnot(gradient_error < 1e-6,
      identical(pc_starts(y, z_PC, M), pc_starts(y, z_altered, M)),
      identical(pc_value_gradient(initial, y, z_PC, M),
                pc_value_gradient(initial, y, z_altered, M)))
    cat('Analytic-gradient check:', gradient_error,
        '\nMasked-label invariance check: PASS\n')
  }
  pc <- fit_PC(y = y, z = z_PC, M = M)

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
        convergence_PC = NA_integer_,
        pc_max_abs_gradient = NA_real_, pc_min_hessian_eigenvalue = NA_real_,
        pc_reproduced_starts = NA_integer_
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
      convergence_PC = pc$convergence,
      pc_max_abs_gradient = pc$max_gradient,
      pc_min_hessian_eigenvalue = pc$min_eigen,
      pc_reproduced_starts = pc$reproduced_starts
    )
  )
}

# ======================================================================
# 14. Versioned checkpointing: old results cannot be reused
config <- list(version = run_version, n = n, B = B, seed = base_seed,
               pi1 = pi1_true, Delta = Delta_true, gamma = gamma_target, a = a_true,
               rng = RNGkind())
checkpoint_file <- file.path(out_dir, 'unequal_prior_OBSERVED_LABEL_ONLY_checkpoint.rds')
sim_results <- data.frame()
if (file.exists(checkpoint_file)) {
  saved <- readRDS(checkpoint_file)
  if (!is.list(saved) || !identical(saved$config, config)) {
    stop('Checkpoint configuration differs. Use a fresh output folder.')
  }
  sim_results <- saved$results
  if (nrow(sim_results)) {
    stopifnot(!anyDuplicated(sim_results$rep),
              all(sim_results$rep %in% seq_len(B)), !anyNA(sim_results$success))
  }
}
completed <- if (nrow(sim_results)) sim_results$rep[sim_results$success] else integer(0)
remaining <- setdiff(seq_len(B), completed)
cat('\nValidated replications already saved:', length(completed),
    '\nReplications to run:', length(remaining), '\n')
for (b in remaining) {
  result <- simulate_one(b)
  if (nrow(sim_results)) sim_results <- sim_results[sim_results$rep != b, , drop = FALSE]
  sim_results <- rbind(sim_results, result)
  if (b %% 25L == 0L || b == B || !result$success) {
    cat(sprintf('Replication %d / %d | validated = %d\n',
                b, B, sum(sim_results$success)))
    flush.console()
  }
  if (b %% 50L == 0L || b == B || !result$success) {
    # Replace only a checkpoint from THIS corrected version.
    saveRDS(list(config = config, results = sim_results), checkpoint_file)
  }
  if (!isTRUE(result$success)) {
    write_csv(sim_results, file.path(out_dir, 'unequal_prior_fit_failure.csv'))
    stop('A replication failed validation. Results have been saved; send the failure CSV. ',
         'No risk summary conditional on successful fits has been produced.')
  }
}
sim_results <- sim_results[order(sim_results$rep), , drop = FALSE]
stopifnot(nrow(sim_results) == B, identical(as.integer(sim_results$rep), seq_len(B)),
          all(sim_results$success), all(sim_results$convergence_CC == 0L),
          all(sim_results$convergence_PC == 0L),
          all(sim_results$pc_max_abs_gradient <= 1e-7),
          all(sim_results$pc_min_hessian_eigenvalue > 1e-8),
          all(sim_results$pc_reproduced_starts >= 2L))

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
  
  geom_ribbon(
    data = risk_summary,
    aes(x = log_c, ymin = n * lower_95, ymax = n * upper_95),
    inherit.aes = FALSE, fill = "#0072B2", alpha = 0.16
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
    title = "n = 2,000; 3,000 paired replications",
    caption = "Blue ribbon: 95% pointwise Monte Carlo intervals from paired risk differences.",
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
    plot.title = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 0),
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
cat("\nSEND THESE THREE FILES FOR THE FINAL MANUSCRIPT CHECK:\n",
    "1. unequal_prior_raw_results_CORRECTED.csv\n",
    "2. unequal_prior_risk_summary_CORRECTED.csv\n",
    "3. unequal_prior_final_summary_CORRECTED.csv\n", sep = "")
