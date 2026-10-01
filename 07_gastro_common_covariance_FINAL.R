# Gastrointestinal data: common-covariance CC versus informative FULL.
# Base R only. No new data, penalty, covariance floor or label imputation.
# NEW specification: two unrestricted Gaussian means with ONE shared covariance.
# FULL missingness: q=expit(xi0+xi1*log(binary posterior entropy)).
# Six training-only starts; two additional balanced-prior starts if checks fail.
# Two optimizers, analytic gradients and Hessian checks; base R only.
# Features are the four previously selected features, held fixed here; this
# analysis does NOT externally validate the historical feature-selection step.
# A numerical pass means a checked local optimum, not proven global optimality.
# Run: setwd("YOUR PROJECT DIRECTORY"); source("07_gastro_common_covariance_FINAL.R")
# Or: Rscript 07_gastro_common_covariance_FINAL.R PROJECT_DIR INPUT_RDS [CHECKPOINT_RDS]
# Optional third argument imports the already validated common-covariance results.
# Otherwise the existing review_gastro_common_covariance checkpoint is reused.
# Every saved fit is rechecked against the input before its prediction is reused.
# Script 06 and the manual fold-53 refinement are not prerequisites.
# No random initialization is used. Platform-dependent optimizer differences may
# leave a fold unresolved; in that event no complete risk curve is exported.

args <- commandArgs(trailingOnly=TRUE)
project_dir <- if(length(args)) args[1L] else getwd()
input_file <- if(length(args)>=2L) args[2L] else NULL
VERSION <- "common_covariance_final_v2"
GRAD_TOL <- 1e-5; OBJECTIVE_TOL <- 1e-6; HESS_TOL <- 1e-7

read_input <- function(f) {
  o <- readRDS(f)
  list(X=as.matrix(o$X), zm=if(!is.null(o$zm)) o$zm else o$z_pc,
       truth=if(!is.null(o$labeltruth)) o$labeltruth else
         if(!is.null(o$truth)) o$truth else o$z_true)
}
if(is.null(input_file)) {
  files <- list.files(project_dir,
    pattern="^gastro_decision_target_analysis(\\([0-9]+\\))?\\.rds$",
    full.names=TRUE, ignore.case=TRUE)
  if(!length(files)) stop("No gastro_decision_target_analysis*.rds found. Put the supplied RDS in the working directory, or give its full path as the second Rscript argument.")
  files <- files[order(file.info(files)$mtime, decreasing=TRUE)]
  inputs <- lapply(files, read_input)
  if(!all(vapply(inputs, function(x) isTRUE(all.equal(x, inputs[[1L]],
      check.attributes=FALSE)), logical(1))))
    stop("The matching RDS files contain different data. Specify the intended input file explicitly.")
  input_file <- files[1L]
}
if(!file.exists(input_file)) stop("Input file does not exist: ",input_file)
d <- read_input(input_file); X <- d$X; zm <- as.numeric(d$zm)
truth <- as.integer(d$truth)
stopifnot(nrow(X)==76L, ncol(X)==4L, length(zm)==76L,
          length(truth)==76L, all(is.finite(X)), all(truth %in% 1:2),
          all(zm[!is.na(zm)] %in% 1:2))
if(any(zm[!is.na(zm)]!=truth[!is.na(zm)])) stop("Available labels differ from histological truth.")
out_dir <- file.path(project_dir,"review_gastro_common_covariance")
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
cat("Input:",normalizePath(input_file),"\nOutput:",out_dir,"\n")

softplus <- function(x) pmax(x,0)+log1p(exp(-abs(x)))
log_entropy <- function(eta) {
  # For |eta|>35 the relative tail-approximation error is below 1e-14.
  # This avoids posterior saturation without imposing a fitted entropy floor.
  t <- abs(eta); small <- t<=35
  lh <- log1p(t)-t; derivative <- -sign(eta)*t/(1+t)
  h <- softplus(-t[small])+t[small]*plogis(-t[small])
  lh[small] <- log(h)
  derivative[small] <- -eta[small]*exp(-softplus(t[small])-
    softplus(-t[small]))/h
  list(value=lh, derivative=derivative)
}
pack <- function(a,mu,S,xi=NULL) {
  L <- t(chol(S)); idx <- lower.tri(L,diag=TRUE)
  v <- L[idx]; diag_idx <- row(L)[idx]==col(L)[idx]
  v[diag_idx] <- log(v[diag_idx])
  c(a,mu[,1L],mu[,2L],v,xi)
}
unpack <- function(theta,p,full=TRUE) {
  k <- p*(p+1L)/2L; L <- matrix(0,p,p)
  idx <- lower.tri(L,diag=TRUE); diagonal <- row(L)[idx]==col(L)[idx]
  v <- theta[seq.int(2L+2L*p,1L+2L*p+k)]
  v[diagonal] <- exp(v[diagonal]); L[idx] <- v
  list(a=theta[1L],mu=matrix(theta[seq.int(2L,1L+2L*p)],p,2L),
       L=L,S=tcrossprod(L),xi=if(full) tail(theta,2L) else NULL,
       idx=idx,diagonal=diagonal)
}
# Negative mean log likelihood and exact analytic gradient.
make_objective <- function(dat,z,full=TRUE) {
  n <- nrow(dat); p <- ncol(dat); missing <- is.na(z)
  cache <- new.env(parent=emptyenv()); cache$theta <- NULL
  evaluate <- function(theta) {
    if(identical(theta,cache$theta)) return(cache$result)
    answer <- tryCatch({
      if(any(!is.finite(theta))) stop("Nonfinite parameter")
      pars <- unpack(theta,p,full); L <- pars$L
      if(any(!is.finite(L)) || any(diag(L)<=0)) stop("Invalid Cholesky factor")
      A <- chol2inv(t(L))
      e1 <- sweep(dat,2L,pars$mu[,1L],"-")
      e2 <- sweep(dat,2L,pars$mu[,2L],"-")
      v1 <- e1 %*% A; v2 <- e2 %*% A
      common <- -.5*(p*log(2*pi)+2*sum(log(diag(L))))
      l1 <- -softplus(-pars$a)+common-.5*rowSums(e1*v1)
      l2 <- -softplus(pars$a)+common-.5*rowSums(e2*v2)
      eta <- l1-l2; tau <- plogis(eta)
      mix <- pmax(l1,l2)+log1p(exp(-abs(eta)))
      r <- tau; r[!missing] <- as.numeric(z[!missing]==1L)
      ll <- mix; ll[!missing & z==1L] <- l1[!missing & z==1L]
      ll[!missing & z==2L] <- l2[!missing & z==2L]
      g <- rep(0,n); gxi <- NULL
      if(full) {
        entropy <- log_entropy(eta)
        predictor <- pars$xi[1L]+pars$xi[2L]*entropy$value
        q <- plogis(predictor); residual <- as.numeric(missing)-q
        ll <- ll+ifelse(missing,-softplus(-predictor),-softplus(predictor))
        g <- residual*pars$xi[2L]*entropy$derivative
        gxi <- c(sum(residual),sum(residual*entropy$value))
      }
      ga <- sum(r-plogis(pars$a))+sum(g)
      gm1 <- colSums(v1*(r+g)); gm2 <- colSums(v2*(1-r-g))
      GS <- .5*(crossprod(v1,v1*(r+g))+crossprod(v2,v2*(1-r-g))-n*A)
      GL <- (GS+t(GS)) %*% L
      gv <- GL[pars$idx]; gv[pars$diagonal] <- gv[pars$diagonal]*diag(L)
      gradient <- -c(ga,gm1,gm2,gv,gxi)/n
      value <- -mean(ll)
      if(!is.finite(value) || any(!is.finite(gradient))) stop("Nonfinite likelihood")
      list(value=value,gradient=gradient,eta=eta,pars=pars)
    },error=function(e) list(value=Inf,gradient=rep(NA_real_,length(theta))))
    cache$theta <- theta; cache$result <- answer; answer
  }
  list(fn=function(theta) evaluate(theta)$value,
       gr=function(theta) evaluate(theta)$gradient, evaluate=evaluate)
}
complete_mle <- function(dat,z) {
  counts <- tabulate(z,2L)
  if(any(counts==0)) stop("CC training class absent")
  mu <- sapply(1:2,function(k) colMeans(dat[z==k,,drop=FALSE]))
  residual <- dat-t(mu[,z,drop=FALSE])
  S <- crossprod(residual)/nrow(dat); chol(S)
  pack(qlogis(counts[1L]/sum(counts)),mu,S)
}
training_starts <- function(dat,z) {
  observed <- !is.na(z); counts <- tabulate(z[observed],2L)
  if(any(counts==0)) stop("FULL has a training class with no available label")
  mu <- sapply(1:2,function(k) colMeans(dat[observed & z==k,,drop=FALSE]))
  e <- dat[observed,,drop=FALSE]-t(mu[,z[observed],drop=FALSE])
  pooled <- crossprod(e)/sum(observed); feature_cov <- cov(dat)
  chol(feature_cov)
  if(inherits(try(chol(pooled),silent=TRUE),"try-error")) pooled <- feature_cov
  blended <- .8*pooled+.2*feature_cov
  shrunk <- .8*mu+.2*matrix(colMeans(dat),ncol(dat),2L)
  configurations <- list(
    list(mu=mu,S=pooled,a=qlogis(counts[1L]/sum(counts))),
    list(mu=mu,S=feature_cov,a=qlogis(counts[1L]/sum(counts))),
    list(mu=mu,S=blended,a=0),
    list(mu=shrunk,S=blended,a=qlogis(counts[1L]/sum(counts))),
    list(mu=mu+matrix(c(.10,-.10,.05,-.05,-.10,.10,-.05,.05),4L,2L),
         S=feature_cov,a=qlogis(counts[1L]/sum(counts))+.2),
    list(mu=mu+matrix(c(-.10,.10,-.05,.05,.10,-.10,.05,-.05),4L,2L),
         S=blended,a=qlogis(counts[1L]/sum(counts))-.2))
  lapply(configurations,function(s) {
    theta <- pack(s$a,s$mu,s$S,c(0,0))
    lh <- log_entropy(make_objective(dat,z)$evaluate(theta)$eta)$value
    initial <- try(suppressWarnings(glm(as.numeric(!observed) ~ lh,
                                       family=binomial())),silent=TRUE)
    xi <- if(!inherits(initial,"try-error") && all(is.finite(coef(initial))))
      as.numeric(coef(initial)) else c(qlogis(mean(!observed)),.5)
    # Clipping applies ONLY to initialization, never to fitted parameters.
    xi <- pmax(-10,pmin(10,xi)); tail_idx <- (length(theta)-1L):length(theta)
    theta[tail_idx] <- xi; theta
  })
}
check_gradient <- function(obj,theta) {
  h <- 1e-5*pmax(1,abs(theta))
  fd <- vapply(seq_along(theta),function(j) {
    plus <- minus <- theta; plus[j] <- plus[j]+h[j]; minus[j] <- minus[j]-h[j]
    (obj$fn(plus)-obj$fn(minus))/(2*h[j])
  },numeric(1))
  error <- max(abs(fd-obj$gr(theta))/pmax(1,abs(fd)))
  if(!is.finite(error) || error>1e-5) stop("Analytic gradient check failed: ",error)
  error
}
hessian_min <- function(obj,theta,step) {
  h <- try(optimHess(theta,obj$fn,obj$gr,
                    control=list(ndeps=rep(step,length(theta)))),silent=TRUE)
  if(inherits(h,"try-error") || any(!is.finite(h))) return(NA_real_)
  min(eigen((h+t(h))/2,symmetric=TRUE,only.values=TRUE)$values)
}
inspect_fit <- function(obj,theta,code,fold,start,stage,message="") {
  e <- obj$evaluate(theta); grad <- max(abs(e$gradient))
  mineig <- cond <- NA_real_
  if(is.finite(e$value)) {
    ev <- eigen(e$pars$S,symmetric=TRUE,only.values=TRUE)$values
    mineig <- min(ev); cond <- max(ev)/min(ev)
  }
  h1 <- h2 <- NA_real_
  if(is.finite(e$value) && is.finite(grad) && grad<=GRAD_TOL && code==0L) {
    h1 <- hessian_min(obj,theta,1e-4); h2 <- hessian_min(obj,theta,5e-5)
  }
  good <- isTRUE(code==0L && is.finite(grad) && grad<=GRAD_TOL &&
    is.finite(mineig) && mineig>0 && is.finite(h1) && is.finite(h2) &&
    min(h1,h2)>HESS_TOL)
  list(theta=theta,record=data.frame(fold=fold,start=start,stage=stage,
    convergence=code,mean_negative_loglik=e$value,max_abs_gradient=grad,
    hessian_min_step1=h1,hessian_min_step2=h2,covariance_min_eigenvalue=mineig,
    covariance_condition=cond,eligible=good,message=message,stringsAsFactors=FALSE))
}
select_candidates <- function(candidates,starts=NULL) {
  records <- do.call(rbind,lapply(candidates,`[[`,"record"))
  finite <- which(is.finite(records$mean_negative_loglik))
  eligible <- which(records$eligible)
  raw <- if(length(finite)) finite[which.min(records$mean_negative_loglik[finite])] else NA_integer_
  selected <- if(length(eligible)) eligible[which.min(records$mean_negative_loglik[eligible])] else raw
  agree <- 0L; lower <- FALSE
  if(length(eligible)) {
    best <- records$mean_negative_loglik[selected]
    near <- eligible[abs(records$mean_negative_loglik[eligible]-best)<=OBJECTIVE_TOL]
    ids <- unique(records$start[near])
    # Agreement must come from genuinely distinct configurations, not merely
    # two optimizers applied to one configuration.
    if(!is.null(starts)) {
      distinct <- list()
      for(id in ids) {
        s <- starts[[as.character(id)]]
        if(is.null(s)) stop("Missing saved starting configuration")
        if(!length(distinct) || !any(vapply(distinct,function(t)
            max(abs(s-t))<1e-10,logical(1)))) distinct[[length(distinct)+1L]] <- s
      }
      agree <- length(distinct)
    } else agree <- length(ids)
    lower <- any(records$mean_negative_loglik[finite]<best-OBJECTIVE_TOL)
  }
  pass <- length(eligible)>0 && agree>=2L && !lower
  reason <- if(!length(eligible)) "No candidate passed convergence, gradient and Hessian checks" else
    if(lower) "An unresolved candidate has a better likelihood" else
    if(agree<2L) "Fewer than two distinct starts agree at the best checked optimum" else
      "Passed numerical local-optimum checks"
  list(pass=pass,reason=reason,agree=agree,records=records,candidates=candidates,
       starting_vectors=starts,
       theta=if(is.na(selected)) NULL else candidates[[selected]]$theta,
       raw_theta=if(is.na(raw)) NULL else candidates[[raw]]$theta)
}
additional_starts <- function(dat,z,original) {
  obj <- make_objective(dat,z)
  starts <- lapply(original[1:2],function(theta) {
    theta[1L] <- 0
    xi_index <- (length(theta)-1L):length(theta)
    theta[xi_index] <- c(0,0)
    lh <- log_entropy(obj$evaluate(theta)$eta)$value
    missing <- as.numeric(is.na(z))
    initial <- try(suppressWarnings(glm(missing ~ lh,family=binomial())),silent=TRUE)
    xi <- if(!inherits(initial,"try-error") && all(is.finite(coef(initial))))
      as.numeric(coef(initial)) else c(qlogis(mean(missing)),.5)
    theta[xi_index] <- pmax(-10,pmin(10,xi)); theta
  })
  names(starts) <- c("7","8"); starts
}
run_starts <- function(obj,starts,fold,candidates=list(),prefix="") {
  for(id in names(starts)) {
    s <- as.integer(id); theta0 <- starts[[id]]
    f <- tryCatch(nlminb(theta0,obj$fn,obj$gr,
      control=list(iter.max=3000L,eval.max=10000L,rel.tol=1e-11,x.tol=1e-10)),
      error=function(e) list(par=theta0,convergence=999L,message=conditionMessage(e)))
    candidates[[length(candidates)+1L]] <- inspect_fit(obj,f$par,f$convergence,
      fold,s,paste0(prefix,"nlminb"),if(is.null(f$message)) "" else f$message)
    begin <- if(is.finite(obj$fn(f$par))) f$par else theta0
    b <- tryCatch(optim(begin,obj$fn,obj$gr,method="BFGS",
      control=list(maxit=3000L,reltol=1e-12)),
      error=function(e) list(par=begin,convergence=999L,message=conditionMessage(e)))
    candidates[[length(candidates)+1L]] <- inspect_fit(obj,b$par,b$convergence,
      fold,s,paste0(prefix,"BFGS"),if(is.null(b$message)) "" else b$message)
  }
  candidates
}
fit_full <- function(dat,z,fold,previous=NULL) {
  obj <- make_objective(dat,z); original <- training_starts(dat,z)
  names(original) <- as.character(seq_along(original))
  if(fold==1L) {
    errors <- vapply(original[1:3],function(s) check_gradient(obj,s),numeric(1))
    cat("FULL analytic-gradient checks:",format(errors),"\n")
  }
  starts <- original
  if(is.null(previous)) candidates <- run_starts(obj,original,fold) else {
    # Recompute likelihoods, gradients, covariance and Hessian checks. Retain
    # optimizer codes and all candidate vectors; never trust a cached pass flag.
    candidates <- lapply(previous$candidates,function(c) {
      rec <- c$record
      if(nrow(rec)!=1L || rec$fold!=fold || !rec$start %in% 1:8)
        stop("Invalid checkpoint candidate identity")
      inspect_fit(obj,c$theta,rec$convergence,fold,rec$start,rec$stage,rec$message)
    })
    if(!length(candidates)) stop("Checkpoint contains no FULL candidates")
    ids <- unique(vapply(candidates,function(c) as.integer(c$record$start),integer(1)))
    if(any(ids>6L)) starts <- c(original,additional_starts(dat,z,original))
    if(!is.null(previous$starting_vectors)) for(id in names(previous$starting_vectors)) {
      if(!id %in% names(starts) || max(abs(previous$starting_vectors[[id]]-starts[[id]]))>1e-8)
        stop("Checkpoint start differs from the specified training-only configuration")
    }
  }
  chosen <- select_candidates(candidates,starts)
  if(!chosen$pass) {
    extra <- additional_starts(dat,z,original)
    already <- unique(chosen$records$start)
    extra <- extra[!as.integer(names(extra)) %in% already]
    # Do not count duplicate initial configurations as independent starts.
    used <- original
    extra <- Filter(function(s) !any(vapply(used,function(t)
      max(abs(s-t))<1e-10,logical(1))),extra)
    if(length(extra)) {
      cat("Fold",fold,": additional training-only balanced-prior starts\n")
      starts <- c(original,additional_starts(dat,z,original))
      candidates <- run_starts(obj,extra,fold,candidates,"additional_")
      chosen <- select_candidates(candidates,starts)
    }
  }
  chosen
}
predict_eta <- function(theta,newdata) {
  if(is.null(theta)) return(NA_real_)
  p <- ncol(newdata); pars <- unpack(theta,p,length(theta)==1+2*p+p*(p+1)/2+2)
  A <- chol2inv(t(pars$L)); e1 <- sweep(newdata,2L,pars$mu[,1L],"-")
  e2 <- sweep(newdata,2L,pars$mu[,2L],"-")
  as.numeric(pars$a-.5*(rowSums((e1%*%A)*e1)-rowSums((e2%*%A)*e2)))
}

checkpoint <- file.path(out_dir,"common_covariance_results.rds")
config <- list(version=VERSION,gradient_tolerance=GRAD_TOL,objective_tolerance=OBJECTIVE_TOL,
               hessian_tolerance=HESS_TOL,n_starts=6L,additional_starts=2L,
               refinement_rule="Use balanced-prior starts 7/8 only after numerical checks fail")
state <- list(config=config,input=d,folds=vector("list",76L),session=sessionInfo())
resume_file <- if(length(args)>=3L) args[3L] else checkpoint
if(file.exists(resume_file)) {
  old <- readRDS(resume_file)
  compatible <- old$config$version %in% c("common_covariance_v1",VERSION) &&
    identical(old$config$gradient_tolerance,GRAD_TOL) &&
    identical(old$config$objective_tolerance,OBJECTIVE_TOL) &&
    identical(old$config$hessian_tolerance,HESS_TOL) &&
    identical(old$config$n_starts,6L) && length(old$folds)==76L &&
    isTRUE(all.equal(old$input,d,check.attributes=FALSE))
  if(!isTRUE(compatible)) stop("Checkpoint data or numerical settings differ from this analysis")
  state <- old
  state$imported_config <- old$config; state$config <- config
  state$verification_session <- sessionInfo()
  # Cached risks are regenerated from the rechecked parameter vectors below.
  state$risk <- NULL; state$ordering <- NULL
  cat("Rechecking saved fits from:",normalizePath(resume_file),"\n")
} else if(length(args)>=3L) stop("Specified checkpoint does not exist: ",resume_file)
# Analytic CC MLE check at a deliberately perturbed parameter vector.
train <- X[-1L,,drop=FALSE]; center <- colMeans(train); scale <- apply(train,2L,sd)
train <- sweep(sweep(train,2L,center,"-"),2L,scale,"/")
cc_test <- complete_mle(train,truth[-1L]); perturb <- seq_along(cc_test)*.001
cat("CC analytic-gradient check:",check_gradient(make_objective(train,truth[-1L],FALSE),cc_test+perturb),"\n")

save_state <- function() {
  available <- Filter(Negate(is.null),state$folds)
  if(length(available)) {
    write.csv(do.call(rbind,lapply(available,`[[`,"prediction")),
              file.path(out_dir,"predictions.csv"),row.names=FALSE)
    write.csv(do.call(rbind,lapply(available,`[[`,"status")),
              file.path(out_dir,"fold_status.csv"),row.names=FALSE)
    write.csv(do.call(rbind,lapply(available,`[[`,"diagnostics")),
              file.path(out_dir,"fit_diagnostics.csv"),row.names=FALSE)
  }
  temp <- paste0(checkpoint,".tmp"); saveRDS(state,temp)
  if(!file.copy(temp,checkpoint,overwrite=TRUE)) stop("Could not save checkpoint")
  unlink(temp)
}
for(i in seq_len(76L)) {
  previous <- state$folds[[i]]
  cat("Fold",i,"/76 ... "); flush.console()
  train <- X[-i,,drop=FALSE]; center <- colMeans(train); scale <- apply(train,2L,sd)
  if(any(!is.finite(scale)) || any(scale<=0)) stop("Training feature has zero/invalid scale in fold ",i)
  train <- sweep(sweep(train,2L,center,"-"),2L,scale,"/")
  test <- sweep(sweep(X[i,,drop=FALSE],2L,center,"-"),2L,scale,"/")
  if(!is.null(previous) && (!isTRUE(all.equal(previous$center,center)) ||
      !isTRUE(all.equal(previous$scale,scale)))) stop("Checkpoint preprocessing differs in fold ",i)
  cc <- complete_mle(train,truth[-i]); cc_obj <- make_objective(train,truth[-i],FALSE)
  cc_diag <- inspect_fit(cc_obj,cc,0L,i,0L,"CC_exact_MLE")
  full <- tryCatch(fit_full(train,zm[-i],i,if(is.null(previous)) NULL else previous$FULL), error=function(e) {
    if(grepl("Analytic gradient check failed",conditionMessage(e),fixed=TRUE)) stop(e)
    record <- cc_diag$record
    record$start <- NA_integer_; record$stage <- "FULL_initialization_error"
    record$convergence <- 999L; record$mean_negative_loglik <- Inf
    record$max_abs_gradient <- NA_real_; record$hessian_min_step1 <- NA_real_
    record$hessian_min_step2 <- NA_real_; record$covariance_min_eigenvalue <- NA_real_
    record$covariance_condition <- NA_real_; record$eligible <- FALSE
    record$message <- conditionMessage(e)
    list(pass=FALSE,reason=conditionMessage(e),agree=0L,records=record,
         theta=NULL,raw_theta=NULL,candidates=list())
  })
  cc_eta <- predict_eta(cc,test); full_eta <- predict_eta(full$theta,test)
  raw_eta <- predict_eta(full$raw_theta,test)
  pass <- cc_diag$record$eligible && full$pass && is.finite(cc_eta) && is.finite(full_eta)
  status <- data.frame(fold=i,CC_pass=cc_diag$record$eligible,FULL_pass=full$pass,
    paired_pass=pass,agreeing_FULL_starts=full$agree,FULL_reason=full$reason)
  prediction <- data.frame(fold=i,truth=truth[i],label_available=!is.na(zm[i]),
    CC_logodds=cc_eta,FULL_selected_logodds=full_eta,FULL_best_raw_logodds=raw_eta,
    FULL_validated_logodds=if(pass) full_eta else NA_real_,paired_pass=pass)
  state$folds[[i]] <- list(center=center,scale=scale,CC_theta=cc,FULL=full,
    status=status,prediction=prediction,diagnostics=rbind(cc_diag$record,full$records))
  save_state(); cat(if(pass) "PASS" else "UNRESOLVED", "\n")
}
pred <- do.call(rbind,lapply(state$folds,`[[`,"prediction"))
state$all_folds_pass <- all(pred$paired_pass)
cat("\nPaired validated folds:",sum(pred$paired_pass),"/76\n")
if(!state$all_folds_pass) {
  # Remove only obsolete outputs from THIS analysis, never the original results.
  obsolete <- c("risk_curve.csv","exact_intervals.csv","exact_transitions.csv",
                "Figure_Gastro_Common_Covariance_FINAL.pdf","Figure_Gastro_Common_Covariance_FINAL.png",
                "entropy_summary.csv")
  unlink(file.path(out_dir,obsolete)); state$risk <- NULL; state$ordering <- NULL
  save_state()
  cat("No risk curve produced: all 76 paired folds must pass. Send the four diagnostic files below.\n")
} else {
  lc <- pred$CC_logodds; lp <- pred$FULL_validated_logodds
  risk_log <- function(log_cost) {
    cc <- lc>log_cost; pc <- lp>log_cost; cost <- exp(log_cost)
    c(CC=(sum(!cc & truth==1L)+cost*sum(cc & truth==2L))/76,
      FULL=(sum(!pc & truth==1L)+cost*sum(pc & truth==2L))/76)
  }
  coefficients_log <- function(log_cost) {
    cc <- lc>log_cost; pc <- lp>log_cost
    c(A=sum(!pc & truth==1L)-sum(!cc & truth==1L),
      B=sum(pc & truth==2L)-sum(cc & truth==2L))
  }
  # Work with the original log-odds at breakpoints. log(exp(eta)) may round
  # differently from eta, which otherwise changes the reported tie convention.
  bounds_log <- sort(unique(c(-3,c(lc,lp)[abs(c(lc,lp))<3],3)))
  pieces <- list()
  for(j in seq_len(length(bounds_log)-1L)) {
    lo <- bounds_log[j]; hi <- bounds_log[j+1L]
    ab <- coefficients_log((lo+hi)/2)
    root <- if(ab[2L]!=0) -ab[1L]/ab[2L] else NA_real_
    root_log <- if(is.finite(root) && root>0) log(root) else NA_real_
    cuts <- c(lo,if(is.finite(root_log) && root_log>lo && root_log<hi) root_log,hi)
    for(k in seq_len(length(cuts)-1L)) {
      left <- cuts[k]; right <- cuts[k+1L]
      pieces[[length(pieces)+1L]] <- data.frame(left=exp(left),right=exp(right),
        log_left=left,log_right=right,A=unname(ab[1L]),B=unname(ab[2L]),
        sign=sign(ab[1L]+exp((left+right)/2)*ab[2L]),n=76L)
    }
  }
  intervals <- do.call(rbind,pieces)
  transitions <- data.frame(c=numeric(),log_c=numeric(),kind=character(),sign_before=numeric(),
    sign_after=numeric(),left_limit=numeric(),right_limit=numeric(),value_at=numeric())
  if(nrow(intervals)>1L) for(j in seq_len(nrow(intervals)-1L)) {
    a <- intervals[j,]; b <- intervals[j+1L,]; cost <- a$right
    if(a$sign!=b$sign) {
      exact <- coefficients_log(a$log_right)
      transitions[nrow(transitions)+1L,] <- list(cost,a$log_right,
        if(a$A==b$A && a$B==b$B) "interior_affine_zero" else "posterior_odds_jump",
        a$sign,b$sign,(a$A+a$B*cost)/76,(b$A+b$B*cost)/76,
        (exact[1L]+exact[2L]*cost)/76)
    }
  }
  grid_log <- sort(unique(c(seq(-3,3,length.out=601L),0,bounds_log,transitions$log_c)))
  risks <- t(vapply(grid_log,risk_log,numeric(2L)))
  curve <- data.frame(c=exp(grid_log),log_c=grid_log,risk_CC=risks[,1L],risk_FULL=risks[,2L],
                      difference=risks[,2L]-risks[,1L])
  # Collapse consecutive pieces with the same ordering for Table S1.
  group <- cumsum(c(TRUE,diff(intervals$sign)!=0))
  ordering <- do.call(rbind,lapply(split(intervals,group),function(s)
    data.frame(left=s$left[1L],right=tail(s$right,1L),
      lower_risk=if(s$sign[1L]<0) "FULL" else if(s$sign[1L]>0) "CC" else "Tie")))
  write.csv(curve,file.path(out_dir,"risk_curve.csv"),row.names=FALSE)
  write.csv(intervals,file.path(out_dir,"exact_intervals.csv"),row.names=FALSE)
  write.csv(transitions,file.path(out_dir,"exact_transitions.csv"),row.names=FALSE)
  entropy <- exp(log_entropy(lp)$value)
  entropy_summary <- do.call(rbind,lapply(c(TRUE,FALSE),function(a) {
    e <- entropy[pred$label_available==a]
    data.frame(label_available=a,n=length(e),mean_entropy=mean(e),median_entropy=median(e))
  }))
  write.csv(entropy_summary,file.path(out_dir,"entropy_summary.csv"),row.names=FALSE)
  # PCA is a whole-data descriptive projection only; it is never used by a fit.
  pca <- prcomp(X,center=TRUE,scale.=TRUE)
  variance <- 100*pca$sdev^2/sum(pca$sdev^2)
  draw_curve <- function(xrange=c(-3,3),inset=FALSE) {
    use <- intervals$log_right>xrange[1L] & intervals$log_left<xrange[2L]
    ss <- intervals[use,,drop=FALSE]
    ys <- c(0,(ss$A+ss$B*pmax(ss$left,exp(xrange[1L])))/76,
      (ss$A+ss$B*pmin(ss$right,exp(xrange[2L])))/76)
    limits <- range(ys); pad <- max(diff(limits)*.08,.001)
    plot(NA,xlim=xrange,ylim=limits+c(-pad,pad),xaxs="i",
      xlab=if(inset) "" else expression(log(c)~","~c==C[12]/C[21]),
      ylab=if(inset) "" else expression(hat(R)[FULL](c)-hat(R)[CC](c)),
      main=if(inset) "Low-cost detail" else "(c) Deployment risk",cex.main=if(inset) .7 else 1)
    for(j in seq_len(nrow(ss))) {
      s <- ss[j,]; x1 <- max(s$log_left,xrange[1L]); x2 <- min(s$log_right,xrange[2L])
      rect(x1,par("usr")[3L],x2,par("usr")[4L],border=NA,
        col=if(s$sign<0) "#E7F2F8" else if(s$sign>0) "#FBEFE5" else "#F0F0F0")
    }
    abline(h=0,lty=2,col="grey45")
    for(j in seq_len(nrow(ss))) {
      s <- ss[j,]; xx <- seq(max(s$log_left,xrange[1L]),min(s$log_right,xrange[2L]),length.out=20L)
      lines(xx,(s$A+s$B*exp(xx))/76,col="#0072B2",lwd=1.5)
    }
    # Vertical segments explicitly mark discontinuities in the empirical risk.
    if(nrow(ss)>1L) for(j in seq_len(nrow(ss)-1L)) {
      a <- ss[j,]; b <- ss[j+1L,]; xx <- a$log_right
      if(xx>xrange[1L] && xx<xrange[2L])
        segments(xx,(a$A+a$B*a$right)/76,xx,(b$A+b$B*b$left)/76,col="#0072B2",lwd=1)
    }
    if(!inset) {
      abline(v=0,lty=3,col="grey45")
      legend("bottomright",c("FULL favoured","CC favoured","Tie"),
        fill=c("#E7F2F8","#FBEFE5","#F0F0F0"),bty="n",cex=.65)
    }
    box()
  }
  draw <- function() {
    par(mfrow=c(1,3),mar=c(4.4,4.1,2.2,.8),mgp=c(2.5,.7,0),cex=.85)
    colour <- ifelse(truth==1L,"#D55E00","#0072B2")
    symbol <- ifelse(pred$label_available,16L,17L)
    plot(pca$x[,1:2],col=colour,pch=symbol,
      xlab=sprintf("PC1 (%.1f%%)",variance[1L]),ylab=sprintf("PC2 (%.1f%%)",variance[2L]),
      main="(a) Feature projection")
    legend("bottomleft",c("Resection","No resection","Available","Masked"),
      col=c("#D55E00","#0072B2","grey30","grey30"),pch=c(16,16,16,17),bty="n",cex=.65)
    boxplot(entropy ~ factor(pred$label_available,levels=c(TRUE,FALSE),labels=c("Available","Masked")),
      col=c("#E7F2F8","#FBEFE5"),outline=FALSE,ylim=c(0,log(2)),
      xlab="Label availability",ylab="Held-out FULL posterior entropy",main="(b) Prediction entropy")
    # Deterministic offsets only affect display; no random-number state is used.
    offset <- .12*sin(seq_len(76L)*2.399963)
    points(ifelse(pred$label_available,1,2)+offset,entropy,col=colour,pch=16,cex=.55)
    draw_curve()
    panel <- par("fig")
    par(fig=c(panel[1L]+.10*diff(panel[1:2]),panel[1L]+.61*diff(panel[1:2]),
              panel[3L]+.59*diff(panel[3:4]),panel[3L]+.91*diff(panel[3:4])),
        new=TRUE,mar=c(1.7,2,1.4,.3),mgp=c(1,.35,0),cex=.55)
    draw_curve(c(-2.04,-1.20),TRUE)
  }
  pdf(file.path(out_dir,"Figure_Gastro_Common_Covariance_FINAL.pdf"),width=12,height=4.4)
  draw(); dev.off()
  png(file.path(out_dir,"Figure_Gastro_Common_Covariance_FINAL.png"),width=7200,height=2640,res=600)
  draw(); dev.off()
  state$risk <- curve
  state$ordering <- list(intervals=intervals,transitions=transitions,summary=ordering)
  state$entropy_summary <- entropy_summary; state$pca_variance_percent <- variance
  save_state()
  cat("Equal-cost risks:\n"); print(risk_log(0))
  cat("Equal-cost confusion counts:\n")
  print(data.frame(method=c("CC","FULL"),false_negative=c(sum(lc<=0 & truth==1L),sum(lp<=0 & truth==1L)),
    false_positive=c(sum(lc>0 & truth==2L),sum(lp>0 & truth==2L))))
  cat("Open intervals (rounded endpoints do not determine endpoint tie decisions):\n"); print(ordering)
  print(entropy_summary)
  cat("Results are exploratory for these fixed features and this shared-covariance model.\n",
      "No independent-fold confidence interval is claimed; overlapping LOOCV training sets are dependent.\n")
}
capture.output(sessionInfo(),file=file.path(out_dir,"sessionInfo.txt"))
cat("\nCompleted. Results and diagnostics are in",normalizePath(out_dir),"\n")
