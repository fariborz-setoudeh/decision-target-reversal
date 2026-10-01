# GitHub entry point: Rscript 02_symmetric_simulation.R PROJECT_DIR
# Run from the repository root in RStudio, or supply its absolute path.
# ======================================================================
# Decision-Target Reversal under Informative Label Missingness
#
# SCRIPT 02: fresh or resumed 5000-replication symmetric simulation
#
# RUN FROM SCRATCH OR EXTEND EXISTING RESULTS TO 5000 REPLICATIONS
# AND PRODUCE FINAL FINITE-SAMPLE / ASYMPTOTIC COMPARISON
#
# ----------------------------------------------------------------------
# This script:
#
# 1. Starts from scratch or reads compatible saved replicate files.
# 2. Preserves all completed replications.
# 3. Adds only the missing replications up to R = 5000.
# 4. Reconstructs the full Monte Carlo summaries.
# 5. Estimates empirical covariance differences D_hat_n.
# 6. Compares:
#
#       actual finite-sample risk difference
#
#       first-order expansion using empirical D_hat_n
#
#       first-order expansion using population D
#
# 7. Estimates empirical critical cost ratios.
# 8. Estimates critical cost ratios implied by D_hat_n.
# 9. Produces final manuscript-quality diagnostic figures.
#
# ======================================================================


# ======================================================================
# 0. CLEAN SESSION
# ======================================================================

rm(list = ls())
gc()

options(
  digits = 12,
  warn = 1
)

set.seed(20260915)


# ======================================================================
# 1. WORKING DIRECTORIES
# ======================================================================

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) args[1L] else getwd()
setwd(project_dir)

results_dir <- file.path(
  getwd(),
  "decision_target_reversal_results"
)

sim_dir <- file.path(
  results_dir,
  "finite_sample_simulation"
)

dir.create(sim_dir, recursive = TRUE, showWarnings = FALSE)


cat(
  "\nWorking directory:\n",
  getwd(),
  "\n",
  sep = ""
)

cat(
  "\nSimulation directory:\n",
  sim_dir,
  "\n",
  sep = ""
)


# ======================================================================
# 2. LOAD VALIDATED POPULATION RESULTS
# ======================================================================

population_file <- file.path(
  results_dir,
  "script01_population_results.RData"
)

if (!file.exists(population_file)) {
  stop("Script 01 population-results file not found.")
}

load(population_file)


required_population_objects <- c(
  "benchmark",
  "benchmark_unequal",
  "transition_table",
  "phase_map"
)

for (obj in required_population_objects) {
  
  if (!exists(obj)) {
    
    stop(
      paste0(
        "Required population object missing: ",
        obj
      )
    )
  }
}


# ======================================================================
# 3. TRUE CONFIGURATION
# ======================================================================

pi_true <- benchmark$pi

Delta_true <- benchmark$Delta

gamma_true <- benchmark$gamma

a_true <- benchmark$a

beta0_true <- benchmark$beta0

xi0_true <- benchmark$xi0

xi1_true <- benchmark$xi1


cat(
  "\n",
  "============================================================\n",
  "TRUE CONFIGURATION\n",
  "============================================================\n",
  sep = ""
)

cat(
  "pi       = ", pi_true, "\n",
  "Delta    = ", Delta_true, "\n",
  "gamma    = ", gamma_true, "\n",
  "a        = ", a_true, "\n",
  "beta0    = ", beta0_true, "\n",
  "xi0      = ", xi0_true, "\n",
  "xi1      = ", xi1_true, "\n",
  "c_minus  = ", benchmark$c_minus, "\n",
  "c_plus   = ", benchmark$c_plus, "\n",
  sep = ""
)


# ======================================================================
# 4. FINAL MONTE CARLO DESIGN
# ======================================================================

TARGET_REP <- 5000L

sample_sizes <- c(
  500L,
  1000L,
  2000L
)


# Same cost grid as Script 02.

cost_grid_base <- exp(
  seq(
    log(0.08),
    log(12),
    length.out = 81L
  )
)

cost_grid <- sort(
  unique(
    c(
      cost_grid_base,
      benchmark$c_minus,
      1,
      benchmark$c_plus
    )
  )
)


cat(
  "\nTarget number of replications per n = ",
  TARGET_REP,
  "\n",
  sep = ""
)


# ======================================================================
# 5. BASIC UTILITIES
# ======================================================================

logit <- function(p) {
  
  log(
    p /
      (1 - p)
  )
}


safe_qlogis <- function(p,
                        eps = 1e-8) {
  
  p <- pmin(
    pmax(
      p,
      eps
    ),
    1 - eps
  )
  
  qlogis(p)
}


logspace_add <- function(a, b) {
  
  m <- pmax(
    a,
    b
  )
  
  m +
    log(
      exp(a - m) +
        exp(b - m)
    )
}


# ======================================================================
# 6. DATA GENERATION
# ======================================================================

generate_sample <- function(n,
                            pi,
                            Delta,
                            xi0,
                            xi1) {
  
  z <- rbinom(
    n,
    size = 1,
    prob = pi
  )
  
  
  mu <- ifelse(
    z == 1,
    Delta / 2,
    -Delta / 2
  )
  
  
  y <- rnorm(
    n,
    mean = mu,
    sd = 1
  )
  
  
  beta0 <- logit(pi)
  
  eta <- beta0 + Delta * y
  
  
  q <- plogis(
    xi0 +
      xi1 * eta^2
  )
  
  
  m <- rbinom(
    n,
    size = 1,
    prob = q
  )
  
  
  data.frame(
    y = y,
    z = z,
    m = m,
    q_true = q
  )
}


# ======================================================================
# 7. COMPLETE-CLASSIFICATION ESTIMATOR
# ======================================================================

fit_CC <- function(dat) {
  
  z <- dat$z
  
  y <- dat$y
  
  
  pi_hat <- mean(z)
  
  
  if (
    pi_hat <= 0 ||
    pi_hat >= 1
  ) {
    
    return(
      list(
        success = FALSE,
        beta0 = NA_real_,
        Delta = NA_real_,
        message = "Class-probability boundary."
      )
    )
  }
  
  
  beta0_hat <- qlogis(
    pi_hat
  )
  
  
  Delta_hat <- 2 *
    mean(
      (2 * z - 1) * y
    )
  
  
  if (
    !is.finite(Delta_hat) ||
    Delta_hat <= 0
  ) {
    
    return(
      list(
        success = FALSE,
        beta0 = beta0_hat,
        Delta = Delta_hat,
        message = "Nonpositive Delta."
      )
    )
  }
  
  
  list(
    success = TRUE,
    beta0 = beta0_hat,
    Delta = Delta_hat,
    message = ""
  )
}


# ======================================================================
# 8. PARTIAL-CLASSIFICATION OBSERVED LOG-LIKELIHOOD
# ======================================================================

negloglik_PC <- function(theta,
                         dat) {
  
  beta0 <- theta[1]
  
  logDelta <- theta[2]
  
  xi0 <- theta[3]
  
  xi1 <- theta[4]
  
  
  if (!all(is.finite(theta))) {
    return(1e100)
  }
  
  
  Delta <- exp(logDelta)
  
  
  if (
    !is.finite(Delta) ||
    Delta < 1e-5 ||
    Delta > 20
  ) {
    
    return(1e100)
  }
  
  
  pi <- plogis(beta0)
  
  
  if (
    pi <= 1e-12 ||
    pi >= 1 - 1e-12
  ) {
    
    return(1e100)
  }
  
  
  y <- dat$y
  
  z <- dat$z
  
  m <- dat$m
  
  
  # ----------------------------------------------------------
  # Marginal density p_beta(y)
  # ----------------------------------------------------------
  
  log1 <-
    log(pi) +
    dnorm(
      y,
      mean = Delta / 2,
      sd = 1,
      log = TRUE
    )
  
  
  log0 <-
    log1p(-pi) +
    dnorm(
      y,
      mean = -Delta / 2,
      sd = 1,
      log = TRUE
    )
  
  
  log_py <- logspace_add(
    log1,
    log0
  )
  
  
  # ----------------------------------------------------------
  # Conditional class-label likelihood
  # ----------------------------------------------------------
  
  eta <- beta0 + Delta * y
  
  
  tau <- plogis(
    eta
  )
  
  
  observed <- (
    m == 0
  )
  
  
  log_label <- numeric(
    length(y)
  )
  
  
  if (any(observed)) {
    
    log_label[observed] <- dbinom(
      z[observed],
      size = 1,
      prob = tau[observed],
      log = TRUE
    )
  }
  
  
  # ----------------------------------------------------------
  # Missingness likelihood
  # ----------------------------------------------------------
  
  q <- plogis(
    xi0 +
      xi1 * eta^2
  )
  
  
  log_missing <- dbinom(
    m,
    size = 1,
    prob = q,
    log = TRUE
  )
  
  
  ll <- sum(
    log_py +
      log_label +
      log_missing
  )
  
  
  if (!is.finite(ll)) {
    return(1e100)
  }
  
  
  -ll
}


# ======================================================================
# 9. PARTIAL-CLASSIFICATION STARTING VALUES
# ======================================================================

initial_PC <- function(dat) {
  
  observed <- (
    dat$m == 0
  )
  
  
  y_obs <- dat$y[
    observed
  ]
  
  
  z_obs <- dat$z[
    observed
  ]
  
  
  if (
    length(z_obs) >= 2 &&
    length(unique(z_obs)) == 2
  ) {
    
    pi_init <- mean(
      z_obs
    )
    
  } else {
    
    pi_init <- 0.5
  }
  
  
  beta0_init <- safe_qlogis(
    pi_init
  )
  
  
  if (
    length(z_obs) >= 2 &&
    length(unique(z_obs)) == 2
  ) {
    
    Delta_init <- 2 *
      mean(
        (2 * z_obs - 1) *
          y_obs
      )
    
  } else {
    
    Delta_init <- 1
  }
  
  
  if (
    !is.finite(Delta_init) ||
    Delta_init <= 0.20
  ) {
    
    Delta_init <- 1
  }
  
  
  Delta_init <- min(
    max(
      Delta_init,
      0.20
    ),
    5
  )
  
  
  eta_init <-
    beta0_init +
    Delta_init * dat$y
  
  
  glm_fit <- try(
    suppressWarnings(
      glm(
        dat$m ~ I(eta_init^2),
        family = binomial()
      )
    ),
    silent = TRUE
  )
  
  
  if (
    !inherits(
      glm_fit,
      "try-error"
    ) &&
    all(
      is.finite(
        coef(glm_fit)
      )
    )
  ) {
    
    xi0_init <- coef(glm_fit)[1]
    
    xi1_init <- coef(glm_fit)[2]
    
  } else {
    
    xi0_init <- safe_qlogis(
      mean(dat$m)
    )
    
    xi1_init <- -0.1
  }
  
  
  xi0_init <- min(
    max(
      xi0_init,
      -5
    ),
    5
  )
  
  
  xi1_init <- min(
    max(
      xi1_init,
      -3
    ),
    1
  )
  
  
  c(
    beta0_init,
    log(Delta_init),
    xi0_init,
    xi1_init
  )
}


# ======================================================================
# 10. FIT PARTIAL-CLASSIFICATION MODEL
# ======================================================================

fit_PC <- function(dat) {
  
  theta0 <- initial_PC(
    dat
  )
  
  
  miss_prop <- mean(
    dat$m
  )
  
  
  miss_prop <- min(
    max(
      miss_prop,
      1e-5
    ),
    1 - 1e-5
  )
  
  
  starts <- rbind(
    
    theta0,
    
    c(
      theta0[1],
      theta0[2],
      qlogis(miss_prop),
      -0.10
    ),
    
    c(
      theta0[1],
      theta0[2],
      qlogis(miss_prop),
      -0.30
    )
  )
  
  
  fits <- vector(
    "list",
    nrow(starts)
  )
  
  
  for (
    s in seq_len(
      nrow(starts)
    )
  ) {
    
    fits[[s]] <- try(
      optim(
        par = starts[s, ],
        fn = negloglik_PC,
        dat = dat,
        method = "BFGS",
        control = list(
          maxit = 1500,
          reltol = 1e-10
        )
      ),
      silent = TRUE
    )
  }
  
  
  valid <- vapply(
    fits,
    function(x) {
      
      !inherits(
        x,
        "try-error"
      ) &&
        is.list(x) &&
        is.finite(x$value) &&
        all(
          is.finite(
            x$par
          )
        )
    },
    logical(1)
  )
  
  
  if (!any(valid)) {
    
    return(
      list(
        success = FALSE,
        beta0 = NA_real_,
        Delta = NA_real_,
        xi0 = NA_real_,
        xi1 = NA_real_,
        convergence = NA_integer_,
        objective = NA_real_,
        start = NA_integer_
      )
    )
  }
  
  
  ids <- which(
    valid
  )
  
  
  values <- vapply(
    fits[ids],
    function(x) x$value,
    numeric(1)
  )
  
  
  best_id <- ids[
    which.min(values)
  ]
  
  
  best <- fits[[
    best_id
  ]]
  
  
  beta0_hat <- best$par[1]
  
  Delta_hat <- exp(
    best$par[2]
  )
  
  xi0_hat <- best$par[3]
  
  xi1_hat <- best$par[4]
  
  
  success <-
    best$convergence == 0 &&
    is.finite(beta0_hat) &&
    is.finite(Delta_hat) &&
    Delta_hat > 0 &&
    is.finite(xi0_hat) &&
    is.finite(xi1_hat)
  
  
  list(
    success = success,
    beta0 = beta0_hat,
    Delta = Delta_hat,
    xi0 = xi0_hat,
    xi1 = xi1_hat,
    convergence = best$convergence,
    objective = best$value,
    start = best_id
  )
}


# ======================================================================
# 11. EXACT POPULATION DEPLOYMENT RISK
# ======================================================================

deployment_risk <- function(beta0_hat,
                            Delta_hat,
                            c,
                            pi_true,
                            Delta_true) {
  
  if (
    !is.finite(beta0_hat) ||
    !is.finite(Delta_hat) ||
    Delta_hat <= 0
  ) {
    
    return(
      rep(
        NA_real_,
        length(c)
      )
    )
  }
  
  
  threshold <-
    (
      log(c) -
        beta0_hat
    ) /
    Delta_hat
  
  
  false_negative <- pnorm(
    threshold,
    mean = Delta_true / 2,
    sd = 1
  )
  
  
  false_positive <- 1 -
    pnorm(
      threshold,
      mean = -Delta_true / 2,
      sd = 1
    )
  
  
  pi_true *
    false_negative +
    c *
    (1 - pi_true) *
    false_positive
}


bayes_risk <- deployment_risk(
  beta0_hat = beta0_true,
  Delta_hat = Delta_true,
  c = cost_grid,
  pi_true = pi_true,
  Delta_true = Delta_true
)


# ======================================================================
# 12. ONE NEW REPLICATION
# ======================================================================

run_one_replication <- function(rep_id,
                                n) {
  
  # EXACTLY SAME SEED RULE AS SCRIPT 02.
  #
  # Thus replications 1:1000 remain unchanged and
  # 1001:5000 form a deterministic continuation.
  
  set.seed(
    20260915 +
      100000L * n +
      rep_id
  )
  
  
  dat <- generate_sample(
    n = n,
    pi = pi_true,
    Delta = Delta_true,
    xi0 = xi0_true,
    xi1 = xi1_true
  )
  
  
  fit_cc <- fit_CC(
    dat
  )
  
  
  fit_pc <- fit_PC(
    dat
  )
  
  
  paired_success <-
    fit_cc$success &&
    fit_pc$success
  
  
  if (paired_success) {
    
    risk_cc <- deployment_risk(
      beta0_hat = fit_cc$beta0,
      Delta_hat = fit_cc$Delta,
      c = cost_grid,
      pi_true = pi_true,
      Delta_true = Delta_true
    )
    
    
    risk_pc <- deployment_risk(
      beta0_hat = fit_pc$beta0,
      Delta_hat = fit_pc$Delta,
      c = cost_grid,
      pi_true = pi_true,
      Delta_true = Delta_true
    )
    
    
    risk_diff <-
      risk_pc -
      risk_cc
    
    
    excess_cc <-
      risk_cc -
      bayes_risk
    
    
    excess_pc <-
      risk_pc -
      bayes_risk
    
  } else {
    
    risk_cc <-
      rep(
        NA_real_,
        length(cost_grid)
      )
    
    
    risk_pc <-
      rep(
        NA_real_,
        length(cost_grid)
      )
    
    
    risk_diff <-
      rep(
        NA_real_,
        length(cost_grid)
      )
    
    
    excess_cc <- risk_cc
    
    excess_pc <- risk_pc
  }
  
  
  parameter_row <- data.frame(
    
    n = n,
    
    rep = rep_id,
    
    paired_success =
      paired_success,
    
    missing_fraction =
      mean(dat$m),
    
    observed_labels =
      sum(dat$m == 0),
    
    cc_success =
      fit_cc$success,
    
    cc_beta0 =
      fit_cc$beta0,
    
    cc_Delta =
      fit_cc$Delta,
    
    pc_success =
      fit_pc$success,
    
    pc_beta0 =
      fit_pc$beta0,
    
    pc_Delta =
      fit_pc$Delta,
    
    pc_xi0 =
      fit_pc$xi0,
    
    pc_xi1 =
      fit_pc$xi1,
    
    pc_convergence =
      fit_pc$convergence,
    
    pc_objective =
      fit_pc$objective,
    
    pc_start =
      fit_pc$start,
    
    stringsAsFactors = FALSE
  )
  
  
  risk_rows <- data.frame(
    
    n = n,
    
    rep = rep_id,
    
    c = cost_grid,
    
    bayes_risk =
      bayes_risk,
    
    risk_CC =
      risk_cc,
    
    risk_PC =
      risk_pc,
    
    excess_CC =
      excess_cc,
    
    excess_PC =
      excess_pc,
    
    diff_PC_minus_CC =
      risk_diff,
    
    paired_success =
      paired_success
  )
  
  
  list(
    parameters = parameter_row,
    risks = risk_rows
  )
}


# ======================================================================
# 13. EXTEND EACH SAMPLE SIZE TO 5000
# ======================================================================

for (
  n_now in sample_sizes
) {
  
  cat(
    "\n",
    "============================================================\n",
    "EXTENDING n = ",
    n_now,
    "\n",
    "============================================================\n",
    sep = ""
  )
  
  
  parameter_file <- file.path(
    sim_dir,
    paste0(
      "replicate_parameters_n",
      n_now,
      ".csv"
    )
  )
  
  
  risk_file <- file.path(
    sim_dir,
    paste0(
      "replicate_risks_n",
      n_now,
      ".csv"
    )
  )
  
  
  if (xor(file.exists(parameter_file), file.exists(risk_file))) {
    stop("Incomplete checkpoint pair for n=", n_now,
         ". Restore both parameter and risk files before resuming.")
  }
  if (file.exists(parameter_file)) {
    parameters_existing <- read.csv(parameter_file)
    risks_existing <- read.csv(risk_file)
    if (anyDuplicated(parameters_existing$rep) ||
        !all(parameters_existing$rep %in% seq_len(TARGET_REP)))
      stop("Invalid or duplicate replication IDs for n=", n_now)
    # A saved parameter row is completed only if its full risk block exists.
    if (!setequal(parameters_existing$rep, unique(risks_existing$rep)))
      stop("Parameter/risk replication IDs differ for n=", n_now)
    for (saved_rep in parameters_existing$rep) {
      cc <- risks_existing$c[risks_existing$rep == saved_rep]
      if (length(cc) != length(cost_grid) || anyDuplicated(cc) ||
          max(abs(sort(cc)-cost_grid)) > 1e-10)
        stop("Incomplete risk grid for n=", n_now, ", rep=", saved_rep)
    }
  } else {
    # Obtain column types without generating or fitting a disposable sample.
    parameters_existing <- data.frame(
      n=integer(),rep=integer(),paired_success=logical(),
      missing_fraction=numeric(),observed_labels=integer(),
      cc_success=logical(),cc_beta0=numeric(),cc_Delta=numeric(),
      pc_success=logical(),pc_beta0=numeric(),pc_Delta=numeric(),
      pc_xi0=numeric(),pc_xi1=numeric(),pc_convergence=integer(),
      pc_objective=numeric(),pc_start=integer())
    risks_existing <- data.frame(n=integer(),rep=integer(),c=numeric(),
      bayes_risk=numeric(),risk_CC=numeric(),risk_PC=numeric(),
      excess_CC=numeric(),excess_PC=numeric(),
      diff_PC_minus_CC=numeric(),paired_success=logical())
  }

  completed <- sort(
    unique(
      parameters_existing$rep
    )
  )
  
  
  cat(
    "Existing completed replications: ",
    length(completed),
    "\n",
    sep = ""
  )
  
  
  reps_to_run <- setdiff(
    seq_len(TARGET_REP),
    completed
  )
  
  
  cat(
    "Additional replications required: ",
    length(reps_to_run),
    "\n",
    sep = ""
  )
  
  
  if (
    length(reps_to_run) > 0
  ) {
    
    buffer_parameters <- list()
    
    buffer_risks <- list()
    
    buffer_count <- 0L
    
    
    for (
      ii in seq_along(
        reps_to_run
      )
    ) {
      
      rep_now <- reps_to_run[ii]
      
      
      ans <- run_one_replication(
        rep_id = rep_now,
        n = n_now
      )
      
      
      buffer_count <- buffer_count + 1L
      
      
      buffer_parameters[[
        buffer_count
      ]] <- ans$parameters
      
      
      buffer_risks[[
        buffer_count
      ]] <- ans$risks
      
      
      # ------------------------------------------------------
      # Save every 50 new replications
      # ------------------------------------------------------
      
      if (
        buffer_count == 50L ||
        ii == length(reps_to_run)
      ) {
        
        add_parameters <- do.call(
          rbind,
          buffer_parameters
        )
        
        
        add_risks <- do.call(
          rbind,
          buffer_risks
        )
        
        
        parameters_existing <- rbind(
          parameters_existing,
          add_parameters
        )
        
        
        risks_existing <- rbind(
          risks_existing,
          add_risks
        )
        
        
        # Remove duplicate replication IDs if the script
        # is restarted after an interrupted save.
        
        parameters_existing <-
          parameters_existing[
            !duplicated(
              parameters_existing$rep,
              fromLast = TRUE
            ),
            ,
            drop = FALSE
          ]
        
        
        risk_key <- paste(
          risks_existing$rep,
          signif(
            risks_existing$c,
            15
          ),
          sep = "_"
        )
        
        
        risks_existing <-
          risks_existing[
            !duplicated(
              risk_key,
              fromLast = TRUE
            ),
            ,
            drop = FALSE
          ]
        
        
        parameters_existing <-
          parameters_existing[
            order(
              parameters_existing$rep
            ),
            ,
            drop = FALSE
          ]
        
        
        risks_existing <-
          risks_existing[
            order(
              risks_existing$rep,
              risks_existing$c
            ),
            ,
            drop = FALSE
          ]
        
        
        write.csv(
          parameters_existing,
          parameter_file,
          row.names = FALSE
        )
        
        
        write.csv(
          risks_existing,
          risk_file,
          row.names = FALSE
        )
        
        
        saveRDS(
          list(
            parameters =
              parameters_existing,
            risks =
              risks_existing
          ),
          file = file.path(
            sim_dir,
            paste0(
              "checkpoint_extended_n",
              n_now,
              ".rds"
            )
          )
        )
        
        
        cat(
          "\r",
          "n=",
          n_now,
          ": total completed = ",
          nrow(
            parameters_existing
          ),
          "/",
          TARGET_REP,
          sep = ""
        )
        
        
        buffer_parameters <- list()
        
        buffer_risks <- list()
        
        buffer_count <- 0L
      }
    }
    
    
    cat("\n")
  }
  
  
  final_parameters <- read.csv(
    parameter_file
  )
  
  
  if (
    length(
      unique(
        final_parameters$rep
      )
    ) != TARGET_REP
  ) {
    
    stop(
      paste0(
        "n=",
        n_now,
        " does not contain exactly ",
        TARGET_REP,
        " replications."
      )
    )
  }
  
  
  cat(
    "Final number of replications for n=",
    n_now,
    ": ",
    nrow(final_parameters),
    "\n",
    sep = ""
  )
}


# ======================================================================
# 14. READ ALL 5000-REPLICATION RESULTS
# ======================================================================

parameter_list <- list()

risk_list <- list()


for (
  n_now in sample_sizes
) {
  
  parameter_list[[
    as.character(
      n_now
    )
  ]] <- read.csv(
    file.path(
      sim_dir,
      paste0(
        "replicate_parameters_n",
        n_now,
        ".csv"
      )
    )
  )
  
  
  risk_list[[
    as.character(
      n_now
    )
  ]] <- read.csv(
    file.path(
      sim_dir,
      paste0(
        "replicate_risks_n",
        n_now,
        ".csv"
      )
    )
  )
}


parameters_all <- do.call(
  rbind,
  parameter_list
)


risks_all <- do.call(
  rbind,
  risk_list
)


rownames(parameters_all) <- NULL

rownames(risks_all) <- NULL


write.csv(
  parameters_all,
  file.path(
    sim_dir,
    "all_replicate_parameters_5000.csv"
  ),
  row.names = FALSE
)


write.csv(
  risks_all,
  file.path(
    sim_dir,
    "all_replicate_risks_5000.csv"
  ),
  row.names = FALSE
)


# ======================================================================
# 15. BASIC MONTE CARLO DIAGNOSTICS
# ======================================================================

diagnostic_rows <- list()


for (
  j in seq_along(
    sample_sizes
  )
) {
  
  n_now <- sample_sizes[j]
  
  
  x <- parameters_all[
    parameters_all$n == n_now,
    ,
    drop = FALSE
  ]
  
  
  good <- x$paired_success
  
  
  diagnostic_rows[[j]] <- data.frame(
    
    n =
      n_now,
    
    replications =
      nrow(x),
    
    successful_pairs =
      sum(good),
    
    success_rate =
      mean(good),
    
    mean_missing_fraction =
      mean(
        x$missing_fraction
      ),
    
    sd_missing_fraction =
      sd(
        x$missing_fraction
      ),
    
    mean_CC_beta0 =
      mean(
        x$cc_beta0[good]
      ),
    
    mean_CC_Delta =
      mean(
        x$cc_Delta[good]
      ),
    
    mean_PC_beta0 =
      mean(
        x$pc_beta0[good]
      ),
    
    mean_PC_Delta =
      mean(
        x$pc_Delta[good]
      ),
    
    mean_PC_xi0 =
      mean(
        x$pc_xi0[good]
      ),
    
    mean_PC_xi1 =
      mean(
        x$pc_xi1[good]
      )
  )
}


diagnostic_table <- do.call(
  rbind,
  diagnostic_rows
)


cat(
  "\n",
  "============================================================\n",
  "FINAL SIMULATION DIAGNOSTICS\n",
  "============================================================\n",
  sep = ""
)


print(
  diagnostic_table
)


write.csv(
  diagnostic_table,
  file.path(
    sim_dir,
    "simulation_diagnostics_5000.csv"
  ),
  row.names = FALSE
)


# ======================================================================
# 16. EMPIRICAL COVARIANCE MATRICES
# ======================================================================

empirical_covariance_list <- list()


for (
  n_now in sample_sizes
) {
  
  x <- parameters_all[
    parameters_all$n == n_now &
      parameters_all$paired_success,
    ,
    drop = FALSE
  ]
  
  
  U_CC <- cbind(
    
    sqrt(n_now) *
      (
        x$cc_beta0 -
          beta0_true
      ),
    
    sqrt(n_now) *
      (
        x$cc_Delta -
          Delta_true
      )
  )
  
  
  U_PC <- cbind(
    
    sqrt(n_now) *
      (
        x$pc_beta0 -
          beta0_true
      ),
    
    sqrt(n_now) *
      (
        x$pc_Delta -
          Delta_true
      )
  )
  
  
  VCC_hat <- cov(
    U_CC
  )
  
  
  VPC_hat <- cov(
    U_PC
  )
  
  
  D_hat <-
    VPC_hat -
    VCC_hat
  
  
  D_hat <- (
    D_hat +
      t(D_hat)
  ) / 2
  
  
  eig_D_hat <- sort(
    eigen(
      D_hat,
      symmetric = TRUE,
      only.values = TRUE
    )$values
  )
  
  
  empirical_covariance_list[[
    as.character(
      n_now
    )
  ]] <- list(
    
    n =
      n_now,
    
    VCC =
      VCC_hat,
    
    VPC =
      VPC_hat,
    
    D =
      D_hat,
    
    eigenvalues =
      eig_D_hat
  )
}


# ======================================================================
# 17. CRITICAL COSTS IMPLIED BY A 2x2 D MATRIX
# ======================================================================

critical_costs_from_D <- function(D,
                                  beta0,
                                  Delta,
                                  tol = 1e-12) {
  
  d00 <- D[1, 1]
  
  d01 <- D[1, 2]
  
  d11 <- D[2, 2]
  
  
  disc <-
    d01^2 -
    d00 * d11
  
  
  if (
    disc <= tol
  ) {
    
    return(
      c(
        c_minus = NA_real_,
        c_plus = NA_real_
      )
    )
  }
  
  
  if (
    abs(d11) <= tol
  ) {
    
    if (
      abs(d01) <= tol
    ) {
      
      return(
        c(
          c_minus = NA_real_,
          c_plus = NA_real_
        )
      )
    }
    
    
    t_root <- -d00 /
      (
        2 * d01
      )
    
    
    c_root <- exp(
      beta0 +
        Delta * t_root
    )
    
    
    return(
      c(
        c_minus = c_root,
        c_plus = NA_real_
      )
    )
  }
  
  
  roots_t <- c(
    
    (
      -d01 -
        sqrt(disc)
    ) /
      d11,
    
    (
      -d01 +
        sqrt(disc)
    ) /
      d11
  )
  
  
  roots_c <- sort(
    exp(
      beta0 +
        Delta * roots_t
    )
  )
  
  
  c(
    c_minus =
      roots_c[1],
    
    c_plus =
      roots_c[2]
  )
}


# ======================================================================
# 18. EMPIRICAL-D CRITICAL COST TABLE
# ======================================================================

covariance_rows <- list()


for (
  j in seq_along(
    sample_sizes
  )
) {
  
  n_now <- sample_sizes[j]
  
  
  obj <- empirical_covariance_list[[
    as.character(
      n_now
    )
  ]]
  
  
  VCC_hat <- obj$VCC
  
  VPC_hat <- obj$VPC
  
  D_hat <- obj$D
  
  
  roots_hat <- critical_costs_from_D(
    D = D_hat,
    beta0 = beta0_true,
    Delta = Delta_true
  )
  
  
  covariance_rows[[j]] <- data.frame(
    
    n =
      n_now,
    
    VCC_00 =
      VCC_hat[1, 1],
    
    VCC_01 =
      VCC_hat[1, 2],
    
    VCC_11 =
      VCC_hat[2, 2],
    
    VPC_00 =
      VPC_hat[1, 1],
    
    VPC_01 =
      VPC_hat[1, 2],
    
    VPC_11 =
      VPC_hat[2, 2],
    
    D00 =
      D_hat[1, 1],
    
    D01 =
      D_hat[1, 2],
    
    D11 =
      D_hat[2, 2],
    
    eig_min =
      obj$eigenvalues[1],
    
    eig_max =
      obj$eigenvalues[2],
    
    empirical_D_c_minus =
      roots_hat["c_minus"],
    
    empirical_D_c_plus =
      roots_hat["c_plus"]
  )
}


covariance_table <- do.call(
  rbind,
  covariance_rows
)


cat(
  "\n",
  "============================================================\n",
  "EMPIRICAL COVARIANCE GEOMETRY\n",
  "============================================================\n",
  sep = ""
)


print(
  covariance_table
)


write.csv(
  covariance_table,
  file.path(
    sim_dir,
    "empirical_covariance_check_5000.csv"
  ),
  row.names = FALSE
)


# ======================================================================
# 19. FIRST-ORDER APPROXIMATION FOR ANY D MATRIX
# ======================================================================

first_order_difference <- function(c,
                                   n,
                                   D,
                                   pi,
                                   Delta,
                                   beta0) {
  
  tC <-
    (
      log(c) -
        beta0
    ) /
    Delta
  
  
  Q <-
    D[1, 1] +
    2 *
    D[1, 2] *
    tC +
    D[2, 2] *
    tC^2
  
  
  mC <-
    pi *
    dnorm(
      tC -
        Delta / 2
    )
  
  
  mC /
    (
      2 *
        Delta *
        n
    ) *
    Q
}


# ======================================================================
# 20. MONTE CARLO RISK SUMMARY
# ======================================================================

summary_rows <- list()

row_id <- 1L


for (
  n_now in sample_sizes
) {
  
  x <- risks_all[
    risks_all$n == n_now &
      risks_all$paired_success,
    ,
    drop = FALSE
  ]
  
  
  D_hat_n <- empirical_covariance_list[[
    as.character(
      n_now
    )
  ]]$D
  
  
  for (
    c_now in cost_grid
  ) {
    
    idx <- abs(
      x$c -
        c_now
    ) <
      1e-10
    
    
    d <- x$diff_PC_minus_CC[
      idx
    ]
    
    
    R <- sum(
      is.finite(d)
    )
    
    
    mean_diff <- mean(
      d,
      na.rm = TRUE
    )
    
    
    sd_diff <- sd(
      d,
      na.rm = TRUE
    )
    
    
    mcse <- sd_diff /
      sqrt(R)
    
    
    pop_approx <- first_order_difference(
      c = c_now,
      n = n_now,
      D = benchmark$D,
      pi = pi_true,
      Delta = Delta_true,
      beta0 = beta0_true
    )
    
    
    empirical_D_approx <- first_order_difference(
      c = c_now,
      n = n_now,
      D = D_hat_n,
      pi = pi_true,
      Delta = Delta_true,
      beta0 = beta0_true
    )
    
    
    summary_rows[[row_id]] <- data.frame(
      
      n =
        n_now,
      
      c =
        c_now,
      
      log_c =
        log(c_now),
      
      successful_reps =
        R,
      
      mean_diff =
        mean_diff,
      
      sd_diff =
        sd_diff,
      
      mcse =
        mcse,
      
      CI_lower =
        mean_diff -
        1.96 * mcse,
      
      CI_upper =
        mean_diff +
        1.96 * mcse,
      
      population_theory =
        pop_approx,
      
      empirical_D_theory =
        empirical_D_approx,
      
      n_mean_diff =
        n_now *
        mean_diff,
      
      n_population_theory =
        n_now *
        pop_approx,
      
      n_empirical_D_theory =
        n_now *
        empirical_D_approx
    )
    
    
    row_id <- row_id + 1L
  }
}


risk_summary <- do.call(
  rbind,
  summary_rows
)


write.csv(
  risk_summary,
  file.path(
    sim_dir,
    "finite_sample_risk_summary_5000.csv"
  ),
  row.names = FALSE
)


# ======================================================================
# 21. EMPIRICAL ZERO CROSSINGS OF ACTUAL RISK CURVE
# ======================================================================

linear_crossings <- function(x,
                             y) {
  
  good <- is.finite(x) &
    is.finite(y)
  
  
  x <- x[good]
  
  y <- y[good]
  
  
  ord <- order(x)
  
  x <- x[ord]
  
  y <- y[ord]
  
  
  idx <- which(
    y[-length(y)] *
      y[-1] <
      0
  )
  
  
  if (
    length(idx) == 0
  ) {
    
    return(
      numeric(0)
    )
  }
  
  
  roots <- numeric(
    length(idx)
  )
  
  
  for (
    j in seq_along(
      idx
    )
  ) {
    
    i <- idx[j]
    
    
    roots[j] <-
      x[i] -
      y[i] *
      (
        x[i + 1] -
          x[i]
      ) /
      (
        y[i + 1] -
          y[i]
      )
  }
  
  
  roots
}


crossing_rows <- list()


for (
  j in seq_along(
    sample_sizes
  )
) {
  
  n_now <- sample_sizes[j]
  
  
  tmp <- risk_summary[
    risk_summary$n == n_now,
    ,
    drop = FALSE
  ]
  
  
  empirical_log_roots <- linear_crossings(
    tmp$log_c,
    tmp$mean_diff
  )
  
  
  empirical_roots <- exp(
    empirical_log_roots
  )
  
  
  D_roots <- critical_costs_from_D(
    D =
      empirical_covariance_list[[
        as.character(
          n_now
        )
      ]]$D,
    beta0 =
      beta0_true,
    Delta =
      Delta_true
  )
  
  
  crossing_rows[[j]] <- data.frame(
    
    n =
      n_now,
    
    risk_curve_c_minus =
      ifelse(
        length(empirical_roots) >= 1,
        empirical_roots[1],
        NA_real_
      ),
    
    risk_curve_c_plus =
      ifelse(
        length(empirical_roots) >= 2,
        empirical_roots[
          length(empirical_roots)
        ],
        NA_real_
      ),
    
    empirical_D_c_minus =
      D_roots["c_minus"],
    
    empirical_D_c_plus =
      D_roots["c_plus"],
    
    population_c_minus =
      benchmark$c_minus,
    
    population_c_plus =
      benchmark$c_plus
  )
}


crossing_table <- do.call(
  rbind,
  crossing_rows
)


cat(
  "\n",
  "============================================================\n",
  "CRITICAL COST COMPARISON\n",
  "============================================================\n",
  sep = ""
)


print(
  crossing_table
)


write.csv(
  crossing_table,
  file.path(
    sim_dir,
    "critical_cost_comparison_5000.csv"
  ),
  row.names = FALSE
)


# ======================================================================
# 22. RMSE OF TWO FIRST-ORDER APPROXIMATIONS
# ======================================================================

rmse_rows <- list()


for (
  j in seq_along(
    sample_sizes
  )
) {
  
  n_now <- sample_sizes[j]
  
  
  tmp <- risk_summary[
    risk_summary$n == n_now,
    ,
    drop = FALSE
  ]
  
  
  rmse_population <- sqrt(
    mean(
      (
        tmp$mean_diff -
          tmp$population_theory
      )^2
    )
  )
  
  
  rmse_empirical_D <- sqrt(
    mean(
      (
        tmp$mean_diff -
          tmp$empirical_D_theory
      )^2
    )
  )
  
  
  rmse_rows[[j]] <- data.frame(
    
    n =
      n_now,
    
    RMSE_population_D =
      rmse_population,
    
    RMSE_empirical_D =
      rmse_empirical_D,
    
    ratio =
      rmse_empirical_D /
      rmse_population
  )
}


rmse_table <- do.call(
  rbind,
  rmse_rows
)


cat(
  "\n",
  "============================================================\n",
  "FIRST-ORDER APPROXIMATION RMSE\n",
  "============================================================\n",
  sep = ""
)


print(
  rmse_table
)


write.csv(
  rmse_table,
  file.path(
    sim_dir,
    "first_order_RMSE_5000.csv"
  ),
  row.names = FALSE
)


# ======================================================================
# 23. EQUAL-COST RESULT
# ======================================================================

equal_cost_table <- risk_summary[
  abs(
    risk_summary$c -
      1
  ) <
    1e-12,
  ,
  drop = FALSE
]


cat(
  "\n",
  "============================================================\n",
  "EQUAL-COST TARGET c=1\n",
  "============================================================\n",
  sep = ""
)


print(
  equal_cost_table[
    ,
    c(
      "n",
      "mean_diff",
      "CI_lower",
      "CI_upper",
      "population_theory",
      "empirical_D_theory"
    )
  ]
)


write.csv(
  equal_cost_table,
  file.path(
    sim_dir,
    "equal_cost_results_5000.csv"
  ),
  row.names = FALSE
)


# ======================================================================
# 24. SELECTED COST VALUES
# ======================================================================

selected_costs <- c(
  0.10,
  benchmark$c_minus,
  0.50,
  1.00,
  2.00,
  benchmark$c_plus,
  8.00
)


selected_rows <- list()

id <- 1L


for (
  n_now in sample_sizes
) {
  
  tmp <- risk_summary[
    risk_summary$n == n_now,
    ,
    drop = FALSE
  ]
  
  
  for (
    cc in selected_costs
  ) {
    
    nearest <- which.min(
      abs(
        tmp$c -
          cc
      )
    )
    
    
    selected_rows[[id]] <-
      tmp[
        nearest,
        ,
        drop = FALSE
      ]
    
    
    id <- id + 1L
  }
}


selected_table <- do.call(
  rbind,
  selected_rows
)


write.csv(
  selected_table,
  file.path(
    sim_dir,
    "selected_cost_results_5000.csv"
  ),
  row.names = FALSE
)


# ======================================================================
# 25. FINAL FIGURE 1
#
# SCALED MONTE CARLO RISK DIFFERENCE
# VERSUS POPULATION FIRST-ORDER THEORY
# ======================================================================

png(
  filename = file.path(
    sim_dir,
    "Figure_final_population_theory_5000.png"
  ),
  width = 1900,
  height = 1450,
  res = 220
)


par(
  mar = c(
    5.3,
    5.5,
    1.5,
    1.5
  )
)


ylim1 <- range(
  c(
    risk_summary$n_mean_diff,
    risk_summary$n_population_theory
  ),
  finite = TRUE
)


plot(
  NA,
  xlim = range(
    risk_summary$log_c
  ),
  ylim = ylim1,
  xlab = expression(
    log(C[12] / C[21])
  ),
  ylab = expression(
    n *
      E[
        R[C](hat(beta)[PC]) -
          R[C](hat(beta)[CC])
      ]
  ),
  main = ""
)


abline(
  h = 0,
  lty = 2
)


abline(
  v = log(
    benchmark$c_minus
  ),
  lty = 3
)


abline(
  v = log(
    benchmark$c_plus
  ),
  lty = 3
)


line_types <- c(
  1,
  2,
  4
)


point_types <- c(
  16,
  17,
  15
)


for (
  j in seq_along(
    sample_sizes
  )
) {
  
  n_now <- sample_sizes[j]
  
  
  tmp <- risk_summary[
    risk_summary$n == n_now,
    ,
    drop = FALSE
  ]
  
  
  lines(
    tmp$log_c,
    tmp$n_mean_diff,
    type = "b",
    lty = line_types[j],
    pch = point_types[j],
    cex = 0.35
  )
}


theory_tmp <- risk_summary[
  risk_summary$n == sample_sizes[1],
  ,
  drop = FALSE
]


lines(
  theory_tmp$log_c,
  theory_tmp$n_population_theory,
  lwd = 2.2
)


legend(
  "topright",
  legend = c(
    paste0(
      "Monte Carlo n=",
      sample_sizes
    ),
    "Population first-order theory"
  ),
  lty = c(
    line_types,
    1
  ),
  pch = c(
    point_types,
    NA
  ),
  lwd = c(
    rep(
      1,
      length(sample_sizes)
    ),
    2.2
  ),
  bty = "n"
)


dev.off()


# ======================================================================
# 26. FINAL FIGURE 2
#
# ACTUAL MONTE CARLO VERSUS FIRST-ORDER APPROXIMATION
# USING EMPIRICAL D_hat_n
# ======================================================================

png(
  filename = file.path(
    sim_dir,
    "Figure_final_empirical_D_5000.png"
  ),
  width = 1900,
  height = 1450,
  res = 220
)


par(
  mar = c(
    5.3,
    5.5,
    1.5,
    1.5
  )
)


ylim2 <- range(
  c(
    risk_summary$n_mean_diff,
    risk_summary$n_empirical_D_theory
  ),
  finite = TRUE
)


plot(
  NA,
  xlim = range(
    risk_summary$log_c
  ),
  ylim = ylim2,
  xlab = expression(
    log(C[12] / C[21])
  ),
  ylab = expression(
    n *
      E[
        R[C](hat(beta)[PC]) -
          R[C](hat(beta)[CC])
      ]
  ),
  main = ""
)


abline(
  h = 0,
  lty = 2
)


abline(
  v = log(
    benchmark$c_minus
  ),
  lty = 3
)


abline(
  v = log(
    benchmark$c_plus
  ),
  lty = 3
)


for (
  j in seq_along(
    sample_sizes
  )
) {
  
  n_now <- sample_sizes[j]
  
  
  tmp <- risk_summary[
    risk_summary$n == n_now,
    ,
    drop = FALSE
  ]
  
  
  # Actual Monte Carlo risk difference
  
  lines(
    tmp$log_c,
    tmp$n_mean_diff,
    lty = line_types[j],
    lwd = 1.8
  )
  
  
  # First-order approximation based on empirical covariance D_hat_n
  
  points(
    tmp$log_c,
    tmp$n_empirical_D_theory,
    pch = point_types[j],
    cex = 0.35
  )
}


legend(
  "topright",
  legend = c(
    paste0(
      "Monte Carlo n=",
      sample_sizes
    ),
    paste0(
      "Empirical-D approximation n=",
      sample_sizes
    )
  ),
  lty = c(
    line_types,
    rep(
      NA,
      length(sample_sizes)
    )
  ),
  pch = c(
    rep(
      NA,
      length(sample_sizes)
    ),
    point_types
  ),
  bty = "n"
)


dev.off()


# ======================================================================
# 27. FINAL FIGURE 3
#
# DIRECT THREE-WAY COMPARISON FOR EACH n
#
# Each panel-like output is saved separately to avoid clutter.
# ======================================================================

for (
  n_now in sample_sizes
) {
  
  tmp <- risk_summary[
    risk_summary$n == n_now,
    ,
    drop = FALSE
  ]
  
  
  png(
    filename = file.path(
      sim_dir,
      paste0(
        "Figure_three_way_n",
        n_now,
        "_5000.png"
      )
    ),
    width = 1800,
    height = 1350,
    res = 220
  )
  
  
  ylim_now <- range(
    c(
      tmp$n_mean_diff,
      tmp$n_empirical_D_theory,
      tmp$n_population_theory
    ),
    finite = TRUE
  )
  
  
  plot(
    tmp$log_c,
    tmp$n_mean_diff,
    type = "l",
    lwd = 2,
    ylim = ylim_now,
    xlab = expression(
      log(C[12] / C[21])
    ),
    ylab = expression(
      n *
        E[
          R[C](hat(beta)[PC]) -
            R[C](hat(beta)[CC])
        ]
    ),
    main = paste0(
      "n = ",
      n_now
    )
  )
  
  
  lines(
    tmp$log_c,
    tmp$n_empirical_D_theory,
    lty = 2,
    lwd = 2
  )
  
  
  lines(
    tmp$log_c,
    tmp$n_population_theory,
    lty = 3,
    lwd = 2
  )
  
  
  abline(
    h = 0,
    lty = 4
  )
  
  
  abline(
    v = log(
      benchmark$c_minus
    ),
    lty = 3
  )
  
  
  abline(
    v = log(
      benchmark$c_plus
    ),
    lty = 3
  )
  
  
  legend(
    "topright",
    legend = c(
      "Monte Carlo",
      expression(
        "First order using " *
          hat(D)[n]
      ),
      "Population first-order theory"
    ),
    lty = c(
      1,
      2,
      3
    ),
    lwd = 2,
    bty = "n"
  )
  
  
  dev.off()
}


# ======================================================================
# 28. SAVE COMPLETE FINAL RESULTS
# ======================================================================

save(
  parameters_all,
  risks_all,
  diagnostic_table,
  empirical_covariance_list,
  covariance_table,
  risk_summary,
  crossing_table,
  rmse_table,
  equal_cost_table,
  selected_table,
  cost_grid,
  sample_sizes,
  TARGET_REP,
  file = file.path(
    sim_dir,
    "script02b_final_5000_results.RData"
  )
)


capture.output(
  sessionInfo(),
  file = file.path(
    sim_dir,
    "sessionInfo_script02b_5000.txt"
  )
)


# ======================================================================
# 29. FINAL AUTOMATIC CHECKS
# ======================================================================

cat(
  "\n",
  "============================================================\n",
  "FINAL AUTOMATIC CHECKS\n",
  "============================================================\n",
  sep = ""
)


mean_missing <- mean(
  parameters_all$missing_fraction
)


success_rate <- mean(
  parameters_all$paired_success
)


cat(
  "\nOverall mean missing fraction = ",
  mean_missing,
  "\n",
  sep = ""
)


cat(
  "Overall paired-fit success rate = ",
  success_rate,
  "\n",
  sep = ""
)


if (
  abs(
    mean_missing -
    gamma_true
  ) >
  0.005
) {
  
  warning(
    "Overall missing fraction differs unexpectedly from gamma."
  )
}


if (
  success_rate <
  0.98
) {
  
  warning(
    "Paired fitting success below 98%."
  )
}


cat(
  "\nPopulation D:\n"
)


print(
  benchmark$D
)


cat(
  "\nEmpirical D matrices:\n"
)


for (
  n_now in sample_sizes
) {
  
  cat(
    "\nn = ",
    n_now,
    "\n",
    sep = ""
  )
  
  
  print(
    empirical_covariance_list[[
      as.character(
        n_now
      )
    ]]$D
  )
}


cat(
  "\nCritical-cost comparison:\n"
)


print(
  crossing_table
)


cat(
  "\nRMSE comparison:\n"
)


print(
  rmse_table
)


cat(
  "\nEqual-cost comparison:\n"
)


print(
  equal_cost_table[
    ,
    c(
      "n",
      "mean_diff",
      "CI_lower",
      "CI_upper",
      "population_theory",
      "empirical_D_theory"
    )
  ]
)


# ======================================================================
# 30. COMPLETION
# ======================================================================

cat(
  "\n",
  "============================================================\n",
  "SCRIPT 02b COMPLETED SUCCESSFULLY\n",
  "============================================================\n",
  sep = ""
)


cat(
  "\nFINAL FILES TO SEND ME:\n\n"
)


cat(
  "1. simulation_diagnostics_5000.csv\n",
  "2. empirical_covariance_check_5000.csv\n",
  "3. finite_sample_risk_summary_5000.csv\n",
  "4. critical_cost_comparison_5000.csv\n",
  "5. first_order_RMSE_5000.csv\n",
  "6. equal_cost_results_5000.csv\n",
  "7. selected_cost_results_5000.csv\n",
  "8. Figure_final_population_theory_5000.png\n",
  "9. Figure_final_empirical_D_5000.png\n",
  "10. Figure_three_way_n500_5000.png\n",
  "11. Figure_three_way_n1000_5000.png\n",
  "12. Figure_three_way_n2000_5000.png\n",
  "13. script02b_final_5000_results.RData\n",
  sep = ""
)


cat(
  "\n============================================================\n"
)