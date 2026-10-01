# GitHub entry point: Rscript 01_population_theory.R PROJECT_DIR
# Run from the repository root in RStudio, or supply its absolute path.
# ======================================================================
# Decision-Target Reversal under Informative Label Missingness
#
# SCRIPT 01 -- 
#
# Population information, geometry transition, and critical cost frontier
#
# ----------------------------------------------------------------------
# MODEL
#
#   Z in {0,1}
#   P(Z=1) = pi
#
#   Y | Z=1 ~ N( Delta/2, 1 )
#   Y | Z=0 ~ N(-Delta/2, 1 )
#
# Posterior log odds:
#
#   eta(y) = beta0 + Delta*y,
#   beta0  = log{pi/(1-pi)}.
#
# Missingness:
#
#   M = 1 means CLASS LABEL IS MISSING
#
#   q(y)
#   = P(M=1 | Y=y)
#   = expit{ xi0 - a * eta(y)^2 }
#
# where xi1 = -a < 0.
#
# For each true configuration, xi0 is calibrated so that
#
#   E[q(Y)] = gamma.
#
# IMPORTANT:
# This calibration is ONLY used to compare data-generating mechanisms
# with the same marginal missing-label proportion gamma.
#
# It is NOT imposed as a likelihood constraint.
# xi0 and xi1 remain nuisance parameters in the statistical model.
#
# ======================================================================


# ======================================================================
# 0. CLEAN SESSION AND WORKING DIRECTORY
# ======================================================================

rm(list = ls())
gc()

options(
  digits = 12,
  warn = 1
)

set.seed(20260914)


# ------------------------------------------------------------
# Set project working directory
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) args[1L] else getwd()
setwd(project_dir)

cat("\nWorking directory:\n")
print(getwd())


# ------------------------------------------------------------
# Output directory
# ------------------------------------------------------------

out_dir <- file.path(
  getwd(),
  "decision_target_reversal_results"
)

if (!dir.exists(out_dir)) {
  dir.create(
    out_dir,
    recursive = TRUE
  )
}

cat("\nOutput directory:\n")
print(out_dir)


# ======================================================================
# 1. BASIC FUNCTIONS
# ======================================================================

logit <- function(p) {
  log(p / (1 - p))
}


# ------------------------------------------------------------
# Marginal density of Y
# ------------------------------------------------------------

fy <- function(y, pi, Delta) {
  
  pi * dnorm(
    y,
    mean = Delta / 2,
    sd = 1
  ) +
    (1 - pi) * dnorm(
      y,
      mean = -Delta / 2,
      sd = 1
    )
}


# ------------------------------------------------------------
# Posterior probability tau(y)=P(Z=1|Y=y)
# ------------------------------------------------------------

tau_y <- function(y, pi, Delta) {
  
  beta0 <- logit(pi)
  
  plogis(
    beta0 + Delta * y
  )
}


# ------------------------------------------------------------
# Posterior log odds
# ------------------------------------------------------------

eta_y <- function(y, pi, Delta) {
  
  logit(pi) + Delta * y
}


# ======================================================================
# 2. EXPECTATIONS UNDER THE MARGINAL MIXTURE
# ======================================================================

# All functions passed to E_mix() MUST be vectorized in y.
#
# integrate() evaluates y in vector form, which is why constructions
# such as c(1,y) must NOT be used inside the integrand.


E_mix <- function(fun,
                  pi,
                  Delta,
                  rel.tol = 1e-9,
                  abs.tol = 1e-10,
                  subdivisions = 2000L) {
  
  integrand <- function(y) {
    
    value <- fun(y)
    
    value * fy(
      y,
      pi = pi,
      Delta = Delta
    )
  }
  
  ans <- integrate(
    f = integrand,
    lower = -Inf,
    upper = Inf,
    rel.tol = rel.tol,
    abs.tol = abs.tol,
    subdivisions = subdivisions,
    stop.on.error = TRUE
  )
  
  ans$value
}


# ======================================================================
# 3. CALIBRATION OF xi0
# ======================================================================

calibrate_xi0 <- function(pi,
                          Delta,
                          gamma,
                          a,
                          tol = 1e-10) {
  
  stopifnot(
    is.finite(pi),
    is.finite(Delta),
    is.finite(gamma),
    is.finite(a),
    
    pi > 0,
    pi < 1,
    
    Delta > 0,
    
    gamma > 0,
    gamma < 1,
    
    a >= 0
  )
  
  
  target <- function(xi0) {
    
    value <- E_mix(
      fun = function(y) {
        
        eta <- eta_y(
          y,
          pi = pi,
          Delta = Delta
        )
        
        plogis(
          xi0 - a * eta^2
        )
      },
      pi = pi,
      Delta = Delta
    )
    
    value - gamma
  }
  
  
  # ----------------------------------------------------------
  # Find a safe bracket automatically
  # ----------------------------------------------------------
  
  lower <- -10
  upper <-  10
  
  
  while (target(lower) > 0) {
    
    lower <- 2 * lower
    
    if (abs(lower) > 1e6) {
      stop("Failed to bracket xi0 from below.")
    }
  }
  
  
  while (target(upper) < 0) {
    
    upper <- 2 * upper
    
    if (abs(upper) > 1e6) {
      stop("Failed to bracket xi0 from above.")
    }
  }
  
  
  uniroot(
    f = target,
    lower = lower,
    upper = upper,
    tol = tol
  )$root
}


# ======================================================================
# 4. POPULATION INFORMATION MATRICES
# ======================================================================

information_matrices <- function(pi,
                                 Delta,
                                 gamma,
                                 a,
                                 rel.tol = 1e-9) {
  
  stopifnot(
    pi > 0,
    pi < 1,
    Delta > 0,
    gamma > 0,
    gamma < 1,
    a >= 0
  )
  
  
  # ----------------------------------------------------------
  # True model parameters
  # ----------------------------------------------------------
  
  beta0 <- logit(pi)
  
  xi1 <- -a
  
  
  # ----------------------------------------------------------
  # Calibrate xi0
  # ----------------------------------------------------------
  
  xi0 <- calibrate_xi0(
    pi = pi,
    Delta = Delta,
    gamma = gamma,
    a = a
  )
  
  
  # ----------------------------------------------------------
  # Missingness probability q(y)
  # ----------------------------------------------------------
  
  qfun <- function(y) {
    
    eta <- beta0 + Delta * y
    
    plogis(
      xi0 + xi1 * eta^2
    )
  }
  
  
  # ----------------------------------------------------------
  # Posterior class probability
  # ----------------------------------------------------------
  
  tfun <- function(y) {
    
    plogis(
      beta0 + Delta * y
    )
  }
  
  
  # ====================================================================
  # 4.1 COMPLETE-CLASSIFICATION INFORMATION
  # ====================================================================
  
  # Under beta = (beta0, Delta),
  #
  # I_CC =
  #
  #   [ pi(1-pi)    0   ]
  #   [    0        1/4  ]
  
  
  ICC <- matrix(
    c(
      pi * (1 - pi), 0,
      0,              1 / 4
    ),
    nrow = 2,
    ncol = 2,
    byrow = TRUE
  )
  
  
  # ====================================================================
  # 4.2 INFORMATION IN Y ALONE
  # ====================================================================
  
  # Marginal score:
  #
  # S_beta0^Y = tau(Y) - pi
  #
  # S_Delta^Y =
  #     Y{2 tau(Y)-1}/2 - Delta/4
  
  
  IY <- matrix(
    0,
    nrow = 2,
    ncol = 2
  )
  
  
  IY[1, 1] <- E_mix(
    fun = function(y) {
      
      tau <- tfun(y)
      
      (tau - pi)^2
    },
    pi = pi,
    Delta = Delta,
    rel.tol = rel.tol
  )
  
  
  IY[1, 2] <- E_mix(
    fun = function(y) {
      
      tau <- tfun(y)
      
      s0 <- tau - pi
      
      s1 <-
        y * (2 * tau - 1) / 2 -
        Delta / 4
      
      s0 * s1
    },
    pi = pi,
    Delta = Delta,
    rel.tol = rel.tol
  )
  
  
  IY[2, 1] <- IY[1, 2]
  
  
  IY[2, 2] <- E_mix(
    fun = function(y) {
      
      tau <- tfun(y)
      
      s1 <-
        y * (2 * tau - 1) / 2 -
        Delta / 4
      
      s1^2
    },
    pi = pi,
    Delta = Delta,
    rel.tol = rel.tol
  )
  
  
  # ====================================================================
  # 4.3 CONDITIONAL LABEL INFORMATION
  # ====================================================================
  
  # I_{Z|Y}(y)
  #
  # = tau(1-tau)
  #
  #     [ 1   y   ]
  #     [ y   y^2 ]
  
  
  IJ <- matrix(
    0,
    nrow = 2,
    ncol = 2
  )
  
  
  IJ[1, 1] <- E_mix(
    fun = function(y) {
      
      tau <- tfun(y)
      
      tau * (1 - tau)
    },
    pi = pi,
    Delta = Delta,
    rel.tol = rel.tol
  )
  
  
  IJ[1, 2] <- E_mix(
    fun = function(y) {
      
      tau <- tfun(y)
      
      tau * (1 - tau) * y
    },
    pi = pi,
    Delta = Delta,
    rel.tol = rel.tol
  )
  
  
  IJ[2, 1] <- IJ[1, 2]
  
  
  IJ[2, 2] <- E_mix(
    fun = function(y) {
      
      tau <- tfun(y)
      
      tau * (1 - tau) * y^2
    },
    pi = pi,
    Delta = Delta,
    rel.tol = rel.tol
  )
  
  
  # ====================================================================
  # 4.4 OBSERVED-LABEL INFORMATION
  # ====================================================================
  
  # Label information is observed only when M=0:
  #
  # E[(1-q(Y)) I_{Z|Y}(Y)]
  
  
  Ilabel_obs <- matrix(
    0,
    nrow = 2,
    ncol = 2
  )
  
  
  Ilabel_obs[1, 1] <- E_mix(
    fun = function(y) {
      
      q   <- qfun(y)
      tau <- tfun(y)
      
      (1 - q) *
        tau *
        (1 - tau)
    },
    pi = pi,
    Delta = Delta,
    rel.tol = rel.tol
  )
  
  
  Ilabel_obs[1, 2] <- E_mix(
    fun = function(y) {
      
      q   <- qfun(y)
      tau <- tfun(y)
      
      (1 - q) *
        tau *
        (1 - tau) *
        y
    },
    pi = pi,
    Delta = Delta,
    rel.tol = rel.tol
  )
  
  
  Ilabel_obs[2, 1] <- Ilabel_obs[1, 2]
  
  
  Ilabel_obs[2, 2] <- E_mix(
    fun = function(y) {
      
      q   <- qfun(y)
      tau <- tfun(y)
      
      (1 - q) *
        tau *
        (1 - tau) *
        y^2
    },
    pi = pi,
    Delta = Delta,
    rel.tol = rel.tol
  )
  
  
  # ====================================================================
  # 4.5 INFORMATION CONTRIBUTED BY MISSINGNESS INDICATOR M
  # ====================================================================
  
  # Missingness linear predictor:
  #
  # xi0 + xi1 * eta(y)^2
  #
  # eta(y) = beta0 + Delta*y
  #
  #
  # beta-score multipliers:
  #
  # d/d beta0:
  #     db0 = 2 xi1 eta
  #
  # d/d Delta:
  #     db1 = 2 xi1 eta y
  #
  #
  # nuisance-score multipliers:
  #
  # d/d xi0:
  #     dxi0 = 1
  #
  # d/d xi1:
  #     dxi1 = eta^2
  #
  #
  # Conditional Bernoulli information weight:
  #
  #     w = q(1-q)
  #
  # IMPORTANT:
  # Every expression below is explicitly vectorized in y.
  
  
  # --------------------------------------------------------------------
  # beta-beta block
  # --------------------------------------------------------------------
  
  Ibb <- matrix(
    0,
    nrow = 2,
    ncol = 2
  )
  
  
  Ibb[1, 1] <- E_mix(
    fun = function(y) {
      
      eta <- beta0 + Delta * y
      
      q <- qfun(y)
      
      w <- q * (1 - q)
      
      db0 <- 2 * xi1 * eta
      
      w * db0^2
    },
    pi = pi,
    Delta = Delta,
    rel.tol = rel.tol
  )
  
  
  Ibb[1, 2] <- E_mix(
    fun = function(y) {
      
      eta <- beta0 + Delta * y
      
      q <- qfun(y)
      
      w <- q * (1 - q)
      
      db0 <- 2 * xi1 * eta
      
      db1 <- 2 * xi1 * eta * y
      
      w * db0 * db1
    },
    pi = pi,
    Delta = Delta,
    rel.tol = rel.tol
  )
  
  
  Ibb[2, 1] <- Ibb[1, 2]
  
  
  Ibb[2, 2] <- E_mix(
    fun = function(y) {
      
      eta <- beta0 + Delta * y
      
      q <- qfun(y)
      
      w <- q * (1 - q)
      
      db1 <- 2 * xi1 * eta * y
      
      w * db1^2
    },
    pi = pi,
    Delta = Delta,
    rel.tol = rel.tol
  )
  
  
  # --------------------------------------------------------------------
  # beta-xi block
  # --------------------------------------------------------------------
  
  Ibxi <- matrix(
    0,
    nrow = 2,
    ncol = 2
  )
  
  
  # beta0 versus xi0
  Ibxi[1, 1] <- E_mix(
    fun = function(y) {
      
      eta <- beta0 + Delta * y
      
      q <- qfun(y)
      
      w <- q * (1 - q)
      
      db0 <- 2 * xi1 * eta
      
      w * db0
    },
    pi = pi,
    Delta = Delta,
    rel.tol = rel.tol
  )
  
  
  # beta0 versus xi1
  Ibxi[1, 2] <- E_mix(
    fun = function(y) {
      
      eta <- beta0 + Delta * y
      
      q <- qfun(y)
      
      w <- q * (1 - q)
      
      db0 <- 2 * xi1 * eta
      
      dxi1 <- eta^2
      
      w * db0 * dxi1
    },
    pi = pi,
    Delta = Delta,
    rel.tol = rel.tol
  )
  
  
  # Delta versus xi0
  Ibxi[2, 1] <- E_mix(
    fun = function(y) {
      
      eta <- beta0 + Delta * y
      
      q <- qfun(y)
      
      w <- q * (1 - q)
      
      db1 <- 2 * xi1 * eta * y
      
      w * db1
    },
    pi = pi,
    Delta = Delta,
    rel.tol = rel.tol
  )
  
  
  # Delta versus xi1
  Ibxi[2, 2] <- E_mix(
    fun = function(y) {
      
      eta <- beta0 + Delta * y
      
      q <- qfun(y)
      
      w <- q * (1 - q)
      
      db1 <- 2 * xi1 * eta * y
      
      dxi1 <- eta^2
      
      w * db1 * dxi1
    },
    pi = pi,
    Delta = Delta,
    rel.tol = rel.tol
  )
  
  
  # --------------------------------------------------------------------
  # xi-xi block
  # --------------------------------------------------------------------
  
  Ixixi <- matrix(
    0,
    nrow = 2,
    ncol = 2
  )
  
  
  # xi0 versus xi0
  Ixixi[1, 1] <- E_mix(
    fun = function(y) {
      
      q <- qfun(y)
      
      q * (1 - q)
    },
    pi = pi,
    Delta = Delta,
    rel.tol = rel.tol
  )
  
  
  # xi0 versus xi1
  Ixixi[1, 2] <- E_mix(
    fun = function(y) {
      
      eta <- beta0 + Delta * y
      
      q <- qfun(y)
      
      w <- q * (1 - q)
      
      w * eta^2
    },
    pi = pi,
    Delta = Delta,
    rel.tol = rel.tol
  )
  
  
  Ixixi[2, 1] <- Ixixi[1, 2]
  
  
  # xi1 versus xi1
  Ixixi[2, 2] <- E_mix(
    fun = function(y) {
      
      eta <- beta0 + Delta * y
      
      q <- qfun(y)
      
      w <- q * (1 - q)
      
      w * eta^4
    },
    pi = pi,
    Delta = Delta,
    rel.tol = rel.tol
  )
  
  
  # ====================================================================
  # 4.6 NUISANCE-PROJECTED EFFICIENT MISSINGNESS INFORMATION
  # ====================================================================
  
  Ixixi_sym <- (
    Ixixi +
      t(Ixixi)
  ) / 2
  
  
  eig_xi <- eigen(
    Ixixi_sym,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  
  
  if (min(eig_xi) <= 1e-10) {
    
    warning(
      paste0(
        "Nuisance information nearly singular. ",
        "Minimum eigenvalue = ",
        format(
          min(eig_xi),
          scientific = TRUE
        )
      )
    )
  }
  
  
  IM_eff <-
    Ibb -
    Ibxi %*%
    solve(
      Ixixi,
      t(Ibxi)
    )
  
  
  IM_eff <- (
    IM_eff +
      t(IM_eff)
  ) / 2
  
  
  # ====================================================================
  # 4.7 EXACT STRUCTURAL CHECK UNDER SYMMETRY
  # ====================================================================
  
  # If pi=1/2, beta0=0 and
  #
  # S_Delta^M
  #   = (2 xi1 / Delta) S_xi1^M.
  #
  # Therefore after nuisance projection:
  #
  #   I_M,22^eff = 0
  #   I_M,12^eff = 0
  #
  # up to numerical precision.
  
  
  if (abs(pi - 0.5) < 1e-12) {
    
    if (
      abs(IM_eff[2, 2]) > 1e-7 ||
      abs(IM_eff[1, 2]) > 1e-7
    ) {
      
      stop(
        paste0(
          "\nSYMMETRY CHECK FAILED.\n",
          "Efficient missingness information in the slope ",
          "direction must vanish at pi=1/2.\n",
          "IM_eff[1,2] = ",
          format(
            IM_eff[1, 2],
            scientific = TRUE
          ),
          "\nIM_eff[2,2] = ",
          format(
            IM_eff[2, 2],
            scientific = TRUE
          ),
          "\n"
        )
      )
    }
  }
  
  
  # ====================================================================
  # 4.8 DIRECT CLOSED-FORM CHECK OF I_M^eff UNDER SYMMETRY
  # ====================================================================
  
  symmetry_closed_form_error <- NA_real_
  
  
  if (abs(pi - 0.5) < 1e-12) {
    
    IM00_direct <- E_mix(
      fun = function(y) {
        
        q <- qfun(y)
        
        w <- q * (1 - q)
        
        4 *
          xi1^2 *
          Delta^2 *
          w *
          y^2
      },
      pi = pi,
      Delta = Delta,
      rel.tol = rel.tol
    )
    
    
    symmetry_closed_form_error <-
      abs(
        IM_eff[1, 1] -
          IM00_direct
      )
  }
  
  
  # ====================================================================
  # 4.9 EFFICIENT INFORMATION UNDER PARTIAL CLASSIFICATION
  # ====================================================================
  
  IPC <-
    IY +
    Ilabel_obs +
    IM_eff
  
  
  IPC <- (
    IPC +
      t(IPC)
  ) / 2
  
  
  # ====================================================================
  # 4.10 INDEPENDENT INFORMATION-DECOMPOSITION CHECK
  # ====================================================================
  
  # Lost conditional-label information
  #
  # E[q(Y) I_{Z|Y}(Y)]
  
  Ilost <-
    IJ -
    Ilabel_obs
  
  
  # Equivalent expression:
  #
  # I_PC
  # =
  # I_CC
  # - E[q I_Z|Y]
  # + I_M^eff
  
  IPC_identity <-
    ICC -
    Ilost +
    IM_eff
  
  
  # ====================================================================
  # 4.11 ASYMPTOTIC COVARIANCES
  # ====================================================================
  
  eig_IPC <- eigen(
    IPC,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  
  
  if (min(eig_IPC) <= 0) {
    
    stop(
      paste0(
        "I_PC is not positive definite. ",
        "Minimum eigenvalue = ",
        min(eig_IPC)
      )
    )
  }
  
  
  VCC <- solve(ICC)
  
  VPC <- solve(IPC)
  
  
  D <-
    VPC -
    VCC
  
  
  D <- (
    D +
      t(D)
  ) / 2
  
  
  eigD <- sort(
    eigen(
      D,
      symmetric = TRUE,
      only.values = TRUE
    )$values
  )
  
  
  # ====================================================================
  # 4.12 NUMERICAL DIAGNOSTICS
  # ====================================================================
  
  gamma_check <- E_mix(
    fun = qfun,
    pi = pi,
    Delta = Delta,
    rel.tol = rel.tol
  )
  
  
  complete_identity_error <-
    max(
      abs(
        ICC -
          (IY + IJ)
      )
    )
  
  
  pc_identity_error <-
    max(
      abs(
        IPC -
          IPC_identity
      )
    )
  
  
  eig_IM_eff <- eigen(
    IM_eff,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  
  
  min_eigen_IMeff <-
    min(eig_IM_eff)
  
  
  min_eigen_IPC <-
    min(eig_IPC)
  
  
  min_eigen_Ixixi <-
    min(eig_xi)
  
  
  # ====================================================================
  # 4.13 CRITICAL COST RATIOS
  # ====================================================================
  
  # For p=1:
  #
  # Q_D(t)
  # =
  # d00 + 2 d01 t + d11 t^2
  #
  # where
  #
  # t_C
  # =
  # {log(c) - beta0}/Delta,
  #
  # c = C12/C21.
  
  
  d00 <- D[1, 1]
  
  d01 <- D[1, 2]
  
  d11 <- D[2, 2]
  
  
  discriminant <-
    d01^2 -
    d00 * d11
  
  
  t_minus <- NA_real_
  
  t_plus <- NA_real_
  
  c_minus <- NA_real_
  
  c_plus <- NA_real_
  
  
  if (
    discriminant > 0 &&
    abs(d11) > 1e-10
  ) {
    
    root1 <-
      (-d01 - sqrt(discriminant)) /
      d11
    
    
    root2 <-
      (-d01 + sqrt(discriminant)) /
      d11
    
    
    t_roots <- sort(
      c(
        root1,
        root2
      )
    )
    
    
    c_roots <- sort(
      exp(
        beta0 +
          Delta * t_roots
      )
    )
    
    
    t_minus <- t_roots[1]
    
    t_plus <- t_roots[2]
    
    c_minus <- c_roots[1]
    
    c_plus <- c_roots[2]
  }
  
  
  # ====================================================================
  # 4.14 RETURN EVERYTHING
  # ====================================================================
  
  list(
    
    # parameters
    pi = pi,
    Delta = Delta,
    gamma = gamma,
    a = a,
    
    beta0 = beta0,
    
    xi0 = xi0,
    xi1 = xi1,
    
    # information
    ICC = ICC,
    
    IY = IY,
    
    IJ = IJ,
    
    Ilabel_obs = Ilabel_obs,
    
    Ilost = Ilost,
    
    Ibb = Ibb,
    
    Ibxi = Ibxi,
    
    Ixixi = Ixixi,
    
    IM_eff = IM_eff,
    
    IPC = IPC,
    
    # covariance
    VCC = VCC,
    
    VPC = VPC,
    
    D = D,
    
    eigD = eigD,
    
    # critical-cost geometry
    discriminant = discriminant,
    
    t_minus = t_minus,
    
    t_plus = t_plus,
    
    c_minus = c_minus,
    
    c_plus = c_plus,
    
    # diagnostics
    diagnostics = c(
      
      gamma_check =
        gamma_check,
      
      complete_identity_error =
        complete_identity_error,
      
      pc_identity_error =
        pc_identity_error,
      
      symmetry_closed_form_error =
        symmetry_closed_form_error,
      
      min_eigen_IMeff =
        min_eigen_IMeff,
      
      min_eigen_IPC =
        min_eigen_IPC,
      
      min_eigen_Ixixi =
        min_eigen_Ixixi
    )
  )
}


# ======================================================================
# 5. Q_D AS A FUNCTION OF THE COST RATIO
# ======================================================================

Q_cost <- function(c, info) {
  
  stopifnot(
    c > 0
  )
  
  
  beta0 <- logit(
    info$pi
  )
  
  
  tC <-
    (
      log(c) -
        beta0
    ) /
    info$Delta
  
  
  D <- info$D
  
  
  D[1, 1] +
    2 * D[1, 2] * tC +
    D[2, 2] * tC^2
}


# ======================================================================
# 6. BENCHMARK CONFIGURATION
# ======================================================================

cat(
  "\n",
  "============================================================\n",
  "BENCHMARK CALCULATION\n",
  "============================================================\n",
  sep = ""
)


benchmark <- information_matrices(
  pi = 0.50,
  Delta = 1.50,
  gamma = 0.30,
  a = 0.25
)


cat("\nParameters\n")
cat("-----------------------------\n")

cat(
  "pi     =",
  benchmark$pi,
  "\n"
)

cat(
  "Delta  =",
  benchmark$Delta,
  "\n"
)

cat(
  "gamma  =",
  benchmark$gamma,
  "\n"
)

cat(
  "a      =",
  benchmark$a,
  "\n"
)

cat(
  "beta0  =",
  benchmark$beta0,
  "\n"
)

cat(
  "xi0    =",
  benchmark$xi0,
  "\n"
)

cat(
  "xi1    =",
  benchmark$xi1,
  "\n"
)


cat("\nI_CC\n")
print(
  benchmark$ICC
)


cat("\nI_Y\n")
print(
  benchmark$IY
)


cat("\nE[I_{Z|Y}]\n")
print(
  benchmark$IJ
)


cat("\nObserved-label information\n")
print(
  benchmark$Ilabel_obs
)


cat("\nLost-label information\n")
print(
  benchmark$Ilost
)


cat("\nI_M^eff\n")
print(
  benchmark$IM_eff
)


cat("\nI_PC^eff\n")
print(
  benchmark$IPC
)


cat("\nV_CC\n")
print(
  benchmark$VCC
)


cat("\nV_PC\n")
print(
  benchmark$VPC
)


cat("\nD = V_PC - V_CC\n")
print(
  benchmark$D
)


cat("\nEigenvalues of D\n")
print(
  benchmark$eigD
)


cat("\nDiscriminant\n")
print(
  benchmark$discriminant
)


cat("\nCritical t values\n")

cat(
  "t_minus =",
  benchmark$t_minus,
  "\n"
)

cat(
  "t_plus  =",
  benchmark$t_plus,
  "\n"
)


cat("\nCritical cost ratios\n")

cat(
  "c_minus =",
  benchmark$c_minus,
  "\n"
)

cat(
  "c_plus  =",
  benchmark$c_plus,
  "\n"
)


cat("\nDiagnostics\n")
print(
  benchmark$diagnostics
)


# ======================================================================
# 7. BENCHMARK ROOT CHECKS
# ======================================================================

if (
  is.finite(benchmark$c_minus) &&
  is.finite(benchmark$c_plus)
) {
  
  cat(
    "\n",
    "============================================================\n",
    "CRITICAL-COST ROOT CHECK\n",
    "============================================================\n",
    sep = ""
  )
  
  
  Q_cminus <-
    Q_cost(
      benchmark$c_minus,
      benchmark
    )
  
  
  Q_cplus <-
    Q_cost(
      benchmark$c_plus,
      benchmark
    )
  
  
  Q_equal <-
    Q_cost(
      1,
      benchmark
    )
  
  
  cat(
    "\nQ(c_minus) =",
    Q_cminus,
    "\n"
  )
  
  
  cat(
    "Q(c_plus)  =",
    Q_cplus,
    "\n"
  )
  
  
  cat(
    "Q(1)       =",
    Q_equal,
    "\n"
  )
  
  
  if (
    abs(Q_cminus) > 1e-7 ||
    abs(Q_cplus) > 1e-7
  ) {
    
    stop(
      "Critical-cost root check failed."
    )
  }
}


# ======================================================================
# 8. ADDITIONAL SYMMETRY CHECK
# ======================================================================

if (
  abs(benchmark$pi - 0.5) < 1e-12
) {
  
  reciprocal_error <-
    abs(
      benchmark$c_minus *
        benchmark$c_plus -
        1
    )
  
  
  cat(
    "\nReciprocal-root check:\n"
  )
  
  
  cat(
    "c_minus * c_plus - 1 =",
    reciprocal_error,
    "\n"
  )
  
  
  if (
    reciprocal_error > 1e-7
  ) {
    
    stop(
      "Symmetric critical-cost reciprocal check failed."
    )
  }
}


# ======================================================================
# 9. SAVE BENCHMARK RESULTS
# ======================================================================

saveRDS(
  benchmark,
  file = file.path(
    out_dir,
    "benchmark_population_information.rds"
  )
)


benchmark_summary <- data.frame(
  
  pi =
    benchmark$pi,
  
  Delta =
    benchmark$Delta,
  
  gamma =
    benchmark$gamma,
  
  a =
    benchmark$a,
  
  xi0 =
    benchmark$xi0,
  
  xi1 =
    benchmark$xi1,
  
  D00 =
    benchmark$D[1, 1],
  
  D01 =
    benchmark$D[1, 2],
  
  D11 =
    benchmark$D[2, 2],
  
  eigen_min =
    benchmark$eigD[1],
  
  eigen_max =
    benchmark$eigD[2],
  
  discriminant =
    benchmark$discriminant,
  
  c_minus =
    benchmark$c_minus,
  
  c_plus =
    benchmark$c_plus,
  
  gamma_check =
    benchmark$diagnostics[
      "gamma_check"
    ],
  
  complete_identity_error =
    benchmark$diagnostics[
      "complete_identity_error"
    ],
  
  pc_identity_error =
    benchmark$diagnostics[
      "pc_identity_error"
    ],
  
  symmetry_closed_form_error =
    benchmark$diagnostics[
      "symmetry_closed_form_error"
    ]
)


write.csv(
  benchmark_summary,
  file = file.path(
    out_dir,
    "benchmark_population_summary.csv"
  ),
  row.names = FALSE
)


# ======================================================================
# 10. FUNCTION FOR MINIMUM EIGENVALUE OF D(a)
# ======================================================================

min_eig_D <- function(a,
                      pi,
                      Delta,
                      gamma) {
  
  info <- information_matrices(
    pi = pi,
    Delta = Delta,
    gamma = gamma,
    a = a
  )
  
  
  min(
    info$eigD
  )
}


# ======================================================================
# 11. LOCATE FIRST TRANSITION a_crit
# ======================================================================

find_acrit <- function(pi,
                       Delta,
                       gamma,
                       a_lower = 0.005,
                       a_upper = 1.00,
                       grid_length = 80L) {
  
  
  stopifnot(
    a_lower > 0,
    a_upper > a_lower
  )
  
  
  grid <- seq(
    from = a_lower,
    to = a_upper,
    length.out = grid_length
  )
  
  
  vals <- vapply(
    X = grid,
    FUN = function(a_now) {
      
      min_eig_D(
        a = a_now,
        pi = pi,
        Delta = Delta,
        gamma = gamma
      )
    },
    FUN.VALUE = numeric(1)
  )
  
  
  # Locate all adjacent sign changes.
  sign_change <- which(
    vals[-length(vals)] *
      vals[-1] <
      0
  )
  
  
  if (
    length(sign_change) == 0
  ) {
    
    return(
      list(
        
        root = NA_real_,
        
        grid = grid,
        
        values = vals
      )
    )
  }
  
  
  # We want the first transition as a increases.
  j <- sign_change[1]
  
  
  root <- uniroot(
    f = function(a_now) {
      
      min_eig_D(
        a = a_now,
        pi = pi,
        Delta = Delta,
        gamma = gamma
      )
    },
    lower = grid[j],
    upper = grid[j + 1],
    tol = 1e-8
  )$root
  
  
  list(
    
    root = root,
    
    grid = grid,
    
    values = vals
  )
}


# ======================================================================
# 12. REPRESENTATIVE TRANSITION CONFIGURATIONS
# ======================================================================

transition_configs <- data.frame(
  
  pi = c(
    0.50,
    0.40,
    0.30,
    0.30,
    0.30,
    0.30
  ),
  
  Delta = c(
    1.50,
    1.50,
    1.50,
    2.00,
    1.50,
    1.50
  ),
  
  gamma = c(
    0.30,
    0.30,
    0.30,
    0.30,
    0.10,
    0.50
  )
)


transition_results <- vector(
  mode = "list",
  length = nrow(
    transition_configs
  )
)


cat(
  "\n",
  "============================================================\n",
  "TRANSITION CALCULATIONS\n",
  "============================================================\n",
  sep = ""
)


for (
  i in seq_len(
    nrow(
      transition_configs
    )
  )
) {
  
  cfg <-
    transition_configs[
      i,
      ,
      drop = FALSE
    ]
  
  
  cat(
    "\nConfiguration ",
    i,
    " of ",
    nrow(
      transition_configs
    ),
    "\n",
    sep = ""
  )
  
  
  cat(
    "pi = ",
    cfg$pi,
    ", Delta = ",
    cfg$Delta,
    ", gamma = ",
    cfg$gamma,
    "\n",
    sep = ""
  )
  
  
  tr <- find_acrit(
    pi = cfg$pi,
    Delta = cfg$Delta,
    gamma = cfg$gamma,
    a_lower = 0.005,
    a_upper = 1.00,
    grid_length = 80L
  )
  
  
  transition_results[[i]] <- data.frame(
    
    pi = cfg$pi,
    
    Delta = cfg$Delta,
    
    gamma = cfg$gamma,
    
    a_crit = tr$root
  )
  
  
  cat(
    "a_crit =",
    tr$root,
    "\n"
  )
}


transition_table <- do.call(
  rbind,
  transition_results
)


cat(
  "\n",
  "============================================================\n",
  "ESTIMATED TRANSITION POINTS\n",
  "============================================================\n\n",
  sep = ""
)


print(
  transition_table
)


write.csv(
  transition_table,
  file = file.path(
    out_dir,
    "theory_transition_points.csv"
  ),
  row.names = FALSE
)


# ======================================================================
# 13. DETAILED PHASE MAP
# ======================================================================

# Main theoretical configuration:
#
# pi     = 0.50
# Delta  = 1.50
# gamma  = 0.30


a_grid <- seq(
  from = 0.02,
  to = 0.60,
  length.out = 100L
)


phase_rows <- vector(
  mode = "list",
  length = length(
    a_grid
  )
)


cat(
  "\n",
  "============================================================\n",
  "COMPUTING SYMMETRIC PHASE MAP\n",
  "============================================================\n",
  sep = ""
)


for (
  i in seq_along(
    a_grid
  )
) {
  
  a_now <- a_grid[i]
  
  
  cat(
    "\r",
    "Phase map: ",
    i,
    "/",
    length(a_grid),
    sep = ""
  )
  
  
  info <- information_matrices(
    pi = 0.50,
    Delta = 1.50,
    gamma = 0.30,
    a = a_now
  )
  
  
  eig <- info$eigD
  
  
  geometry <- if (
    eig[1] > 1e-8 &&
    eig[2] > 1e-8
  ) {
    
    "positive definite"
    
  } else if (
    eig[1] < -1e-8 &&
    eig[2] > 1e-8
  ) {
    
    "indefinite"
    
  } else if (
    eig[1] < -1e-8 &&
    eig[2] < -1e-8
  ) {
    
    "negative definite"
    
  } else {
    
    "boundary"
  }
  
  
  phase_rows[[i]] <- data.frame(
    
    a = a_now,
    
    xi0 = info$xi0,
    
    lambda_min = eig[1],
    
    lambda_max = eig[2],
    
    D00 =
      info$D[1, 1],
    
    D01 =
      info$D[1, 2],
    
    D11 =
      info$D[2, 2],
    
    determinant =
      det(info$D),
    
    discriminant =
      info$discriminant,
    
    c_minus =
      info$c_minus,
    
    c_plus =
      info$c_plus,
    
    log_c_minus =
      ifelse(
        is.finite(
          info$c_minus
        ),
        log(
          info$c_minus
        ),
        NA_real_
      ),
    
    log_c_plus =
      ifelse(
        is.finite(
          info$c_plus
        ),
        log(
          info$c_plus
        ),
        NA_real_
      ),
    
    geometry =
      geometry,
    
    gamma_check =
      info$diagnostics[
        "gamma_check"
      ],
    
    complete_identity_error =
      info$diagnostics[
        "complete_identity_error"
      ],
    
    pc_identity_error =
      info$diagnostics[
        "pc_identity_error"
      ],
    
    symmetry_closed_form_error =
      info$diagnostics[
        "symmetry_closed_form_error"
      ]
  )
}


cat("\n")


phase_map <- do.call(
  rbind,
  phase_rows
)


write.csv(
  phase_map,
  file = file.path(
    out_dir,
    "theory_phase_map_symmetric.csv"
  ),
  row.names = FALSE
)


# ======================================================================
# 14. CHECK IMMEDIATELY BELOW AND ABOVE a_crit
# ======================================================================

index_sym <- which(
  abs(
    transition_table$pi -
      0.50
  ) <
    1e-12 &
    abs(
      transition_table$Delta -
        1.50
    ) <
    1e-12 &
    abs(
      transition_table$gamma -
        0.30
    ) <
    1e-12
)


if (
  length(index_sym) == 1
) {
  
  acrit_sym <-
    transition_table$a_crit[
      index_sym
    ]
  
  
  if (
    is.finite(
      acrit_sym
    )
  ) {
    
    a_test <- c(
      
      0.95 *
        acrit_sym,
      
      acrit_sym,
      
      1.05 *
        acrit_sym
    )
    
    
    cat(
      "\n",
      "============================================================\n",
      "CHECK AROUND SYMMETRIC a_crit\n",
      "============================================================\n",
      sep = ""
    )
    
    
    for (
      aa in a_test
    ) {
      
      info <- information_matrices(
        pi = 0.50,
        Delta = 1.50,
        gamma = 0.30,
        a = aa
      )
      
      
      cat(
        "\na = ",
        format(
          aa,
          digits = 12
        ),
        "\n",
        sep = ""
      )
      
      
      cat(
        "Eigenvalues(D):\n"
      )
      
      
      print(
        info$eigD
      )
      
      
      cat(
        "D:\n"
      )
      
      
      print(
        info$D
      )
      
      
      cat(
        "Critical costs:\n"
      )
      
      
      cat(
        "c_minus =",
        info$c_minus,
        "\n"
      )
      
      
      cat(
        "c_plus  =",
        info$c_plus,
        "\n"
      )
    }
  }
}


# ======================================================================
# 15. UNEQUAL-PRIOR BENCHMARK CHECK
# ======================================================================

cat(
  "\n",
  "============================================================\n",
  "UNEQUAL-PRIOR CHECK\n",
  "============================================================\n",
  sep = ""
)


benchmark_unequal <- information_matrices(
  pi = 0.30,
  Delta = 1.50,
  gamma = 0.30,
  a = 0.25
)


cat(
  "\nxi0 =",
  benchmark_unequal$xi0,
  "\n"
)


cat(
  "\nD:\n"
)


print(
  benchmark_unequal$D
)


cat(
  "\nEigenvalues(D):\n"
)


print(
  benchmark_unequal$eigD
)


cat(
  "\nCritical costs:\n"
)


cat(
  "c_minus =",
  benchmark_unequal$c_minus,
  "\n"
)


cat(
  "c_plus  =",
  benchmark_unequal$c_plus,
  "\n"
)


cat(
  "\nDiagnostics:\n"
)


print(
  benchmark_unequal$diagnostics
)


saveRDS(
  benchmark_unequal,
  file = file.path(
    out_dir,
    "benchmark_unequal_prior.rds"
  )
)


# ======================================================================
# 16. GLOBAL NUMERICAL QUALITY CHECKS FOR PHASE MAP
# ======================================================================

cat(
  "\n",
  "============================================================\n",
  "GLOBAL NUMERICAL CHECKS\n",
  "============================================================\n",
  sep = ""
)


max_gamma_error <-
  max(
    abs(
      phase_map$gamma_check -
        0.30
    ),
    na.rm = TRUE
  )


max_complete_identity_error <-
  max(
    phase_map$complete_identity_error,
    na.rm = TRUE
  )


max_pc_identity_error <-
  max(
    phase_map$pc_identity_error,
    na.rm = TRUE
  )


max_symmetry_error <-
  max(
    phase_map$symmetry_closed_form_error,
    na.rm = TRUE
  )


cat(
  "\nMaximum |gamma_check - gamma| =",
  max_gamma_error,
  "\n"
)


cat(
  "Maximum complete-information identity error =",
  max_complete_identity_error,
  "\n"
)


cat(
  "Maximum PC-information identity error =",
  max_pc_identity_error,
  "\n"
)


cat(
  "Maximum symmetry closed-form error =",
  max_symmetry_error,
  "\n"
)


# Fail rather than silently proceed if something is numerically wrong.

if (
  max_gamma_error > 1e-7
) {
  
  stop(
    "Marginal missingness calibration failed numerical quality check."
  )
}


if (
  max_complete_identity_error > 1e-7
) {
  
  stop(
    "Complete-information decomposition failed."
  )
}


if (
  max_pc_identity_error > 1e-7
) {
  
  stop(
    "Partial-classification information decomposition failed."
  )
}


if (
  max_symmetry_error > 1e-7
) {
  
  stop(
    "Symmetric efficient-information identity failed."
  )
}


# ======================================================================
# 17. SAVE COMPLETE WORKSPACE OBJECTS
# ======================================================================

save(
  benchmark,
  benchmark_unequal,
  transition_table,
  phase_map,
  file = file.path(
    out_dir,
    "script01_population_results.RData"
  )
)


# ======================================================================
# 18. SAVE SESSION INFORMATION
# ======================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    out_dir,
    "sessionInfo_script01.txt"
  )
)


# ======================================================================
# 19. FINAL OUTPUT
# ======================================================================

cat(
  "\n",
  "============================================================\n",
  "SCRIPT 01 COMPLETED SUCCESSFULLY\n",
  "============================================================\n",
  sep = ""
)


cat(
  "\nResults saved in:\n",
  out_dir,
  "\n",
  sep = ""
)


cat(
  "\nFiles created:\n"
)


cat(
  "  1. benchmark_population_summary.csv\n"
)


cat(
  "  2. benchmark_population_information.rds\n"
)


cat(
  "  3. benchmark_unequal_prior.rds\n"
)


cat(
  "  4. theory_transition_points.csv\n"
)


cat(
  "  5. theory_phase_map_symmetric.csv\n"
)


cat(
  "  6. script01_population_results.RData\n"
)


cat(
  "  7. sessionInfo_script01.txt\n"
)


cat(
  "\n============================================================\n"
)