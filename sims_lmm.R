library(tictoc)
library(lme4)
library(lmmstest)
library(tictoc)
# Settings are (n, T, lambda_1, lambda_2, sigma, gamma),
# where gamma is scale factor for lambda
do_one_sim <- function(seed, settings){

  set.seed(seed)

  # Unpack settings
  n_i <- settings[1]
  n_t <- settings[2]
  sig <- settings[3]
  lam_1 <- settings[4] * settings[6]
  lam_2 <- settings[5] * settings[6]

  # Generate data
  X <- matrix(runif(n_i * n_t, -1, 2) , nrow = n_i, ncol = n_t)
  U1 <- rnorm(n_i, sd = lam_1)
  U2 <- rnorm(n_i, sd = lam_2)
  E <- matrix(rnorm(n_i * n_t, sd = sig), nrow = n_i, ncol = n_t)
  Y <- X + 1 + E
  for(ii in 1:n_i){
    Y[ii, ] <- Y[ii, ] + U1[ii] + U2[ii] * X[ii, ]
  }
  D <- data.frame(c(t(Y)), c(t(X)))
  D$time <- rep(1:n_t, n_i)
  D$unit <- as.factor(rep(1:n_i, each = n_t))
  names(D)[1:2] <- c("y", "x")

  # Fit model with sigma unknown
  tictoc::tic()
  fit <- lme4::lmer(y ~ x + (1|unit) + (0 + x|unit), data = D, REML = FALSE)
  timer_mle <- tictoc::toc(quiet = TRUE)
  VC <- as.data.frame(lme4::VarCorr(fit))
  sig_hat <- VC[3, 5]
  lam_hat <- VC[1:2, 5]
  beta_hat <- fixef(fit)

  # Extract model matrices
  Z_tall <- as.matrix(lme4::getME(fit, "Z"))
  X_tall <- lme4::getME(fit, "X")

  # Fit model with sigma known
  tictoc::tic()
  obj <- function(theta){
    -lmmstest::log_lik(y = D$y,
                           X = X_tall,
                           Z = Z_tall,
                           Beta = theta[1:2],
                           sigma = sig,
                           lambda = theta[3:4],
                           lam_idx = rep(1:2, each = n_i),
                           diffs = 0)[[1]]
  }
  obj_grad <- function(theta){
    -lmmstest:::score(y = D$y,
                      X = X_tall,
                      Z = Z_tall,
                      Beta = theta[1:2],
                      sigma = sig,
                      lambda = theta[3:4],
                      lam_idx = rep(1:2, each = n_i))[c(1:2, 4:5)]
  }

  opt <- optim(par = c(beta_hat, lam_hat),
               fn = obj,
               gr = obj_grad,
               method = "L-BFGS-B",
               lower = c(-Inf, -Inf, 0, 0),
               upper = c(Inf, Inf, Inf, Inf),
               control = list(factr = 1e4))
  timer_mle_known <- tictoc::toc(quiet = TRUE)
  #############################################################################
  # For power figure: check if reject (lam_1, lam_2) = (1e-6, 1e-6)
  # NB: sigma is fixed and known
  #############################################################################
  lam_null <- c(1e-6, 1e-6)
  L_null <- diag(rep(lam_null, each = n_i), n_i * 2)
  tictoc::tic()
  Sigma_null <- diag(sig^2, n_i * n_t) + Z_tall %*% L_null^2 %*% t(Z_tall)
  beta_null <- c(qr.solve(crossprod(X_tall, qr.solve(Sigma_null, X_tall)),
                          crossprod(X_tall, qr.solve(Sigma_null, D$y))))
  timer_mle_null <- tictoc::toc(quiet = TRUE)

  # Likelihood ratio test (power)
  tictoc::tic()
  ll_null <- lmmstest::log_lik(y = D$y,
                               X = X_tall,
                               Z = Z_tall,
                               Beta = beta_null,
                               sigma = sig,
                               lambda = lam_null,
                               lam_idx = rep(1:2, each = n_i),
                               diffs = 0)[[1]]
  ll_alt <- lmmstest::log_lik(y = D$y,
                              X = X_tall,
                              Z = Z_tall,
                              Beta = opt$par[1:2],
                              sigma = sig,
                              lambda = opt$par[3:4],
                              lam_idx = rep(1:2, each = n_i),
                              diffs = 0)[[1]]
  lrt_stat_pow <- 2 * (ll_alt - ll_null)
  timer_lrt <- tictoc::toc(quiet = TRUE)
  lrt_p_val_pow <- pchisq(lrt_stat_pow, df = 2, lower.tail = FALSE)

  # Wald test (power)
  tictoc::tic()
  finf <- lmmstest:::fish_inf(y = D$y,
                              X = X_tall,
                              Z = Z_tall,
                              Beta = opt$par[1:2],
                              sigma = sig,
                              lambda = opt$par[3:4],
                              lam_idx = rep(1:2, each = n_i))
  e <- opt$par[3:4] - lam_null
  wald_stat_pow <- c(crossprod(e, finf[4:5, 4:5] %*% e))
  timer_wald <- tictoc::toc(quiet = TRUE)
  wald_p_val_pow <- pchisq(wald_stat_pow, df = 2, lower.tail = FALSE)

  # Proposed test (power)
  tictoc::tic()
  prop_pow <- lmmstest::score_test(
    y = D$y,
    X = X_tall,
    Z = Z_tall,
    Beta = beta_null,
    sigma = sig,
    lambda = lam_null,
    lam_idx = rep(1:2, each = n_i),
    test_idx = 4:5, # lambda parameters
    fix_idx = 3 # Sigma is known
  )
  timer_our <- tictoc::toc(quiet = TRUE)
  #############################################################################
  # For cover figure: (i) check if cover (lam_1, lam_2) with sigma known
  # (ii) check if cover (lam_1) with sigma and lam_2 as nuisance parameters
  #############################################################################
  lam_null <- c(lam_1, lam_2)
  L_null <- diag(rep(lam_null, each = n_i), n_i * 2)
  Sigma_null <- diag(sig^2, n_i * n_t) + Z_tall %*% L_null^2 %*% t(Z_tall)
  beta_null <- c(qr.solve(crossprod(X_tall, qr.solve(Sigma_null, X_tall)),
                          crossprod(X_tall, qr.solve(Sigma_null, D$y))))
  # Likelihood ratio test (cover both)
  ll_null <- lmmstest::log_lik(y = D$y,
                               X = X_tall,
                               Z = Z_tall,
                               Beta = beta_null,
                               sigma = sig,
                               lambda = lam_null,
                               lam_idx = rep(1:2, each = n_i),
                               diffs = 0)[[1]]
  lrt_stat_cov <- 2 * (ll_alt - ll_null)
  lrt_p_val_cov <- pchisq(lrt_stat_cov, df = 2, lower.tail = FALSE)

  # Wald test (cover both)
  e <- opt$par[3:4] - lam_null
  wald_stat_cov<- c(crossprod(e, finf[4:5, 4:5] %*% e))
  wald_p_val_cov <- pchisq(wald_stat_cov, df = 2, lower.tail = FALSE)

  # Proposed test (cover both)
  prop_cov <- lmmstest::score_test(
    y = D$y,
    X = X_tall,
    Z = Z_tall,
    Beta = beta_null,
    sigma = sig,
    lambda = lam_null,
    lam_idx = rep(1:2, each = n_i),
    test_idx = 4:5, # lambda parameters
    fix_idx = 3 # sigma is known
  )

  prop_cov_true <- lmmstest::score_test(
    y = D$y,
    X = X_tall,
    Z = Z_tall,
    Beta = c(1, 1),
    sigma = sig,
    lambda = lam_null,
    lam_idx = rep(1:2, each = n_i),
    test_idx = 4:5, # lambda parameters
    fix_idx = 1:3 # beta and sigma are known
  )

  # Proposed test (cover w sigma nuisance parameter)
  prop_cov_nuis <- lmmstest::score_test(
    y = D$y,
    X = X_tall,
    Z = Z_tall,
    Beta = beta_hat, # Estimate
    sigma = sig_hat, # Estimate
    lambda = c(lam_1, lam_2),
    lam_idx = rep(1:2, each = n_i),
    test_idx = 4:5,
    efficient = TRUE
  )
  # Return
  out <- c(prop_pow, lrt_stat_pow, lrt_p_val_pow,  wald_stat_pow, wald_p_val_pow,
    prop_cov,  lrt_stat_cov, lrt_p_val_cov,  wald_stat_cov, wald_p_val_cov,
    prop_cov_true[-2], prop_cov_nuis[-2], timer_our$toc - timer_our$tic,
    timer_lrt$toc - timer_lrt$tic, timer_wald$toc - timer_wald$tic,
    timer_mle$toc - timer_mle$tic, timer_mle_known$toc - timer_mle_known$tic,
    timer_mle_null$toc - timer_mle_null$tic
    )
    names(out) <- c("stat_our_pow", "df_pow", "pval_our_pow",
      "stat_lrt_pow", "pval_lrt_pow",  "stat_wald_pow", "pval_wald_pow",
      "stat_our", "df", "pval_our",
      "stat_lrt", "pval_lrt", "stat_wald", "pval_wald",
      "stat_our_true", "pval_our_true", "stat_our_nuis", "pval_our_nuis",
      "time_our", "time_lrt", "time_wald", "time_mle", "time_mle_known",
      "time_mle_null")
    return(out)
}

# Do simulation
library(doParallel)
library(doRNG)
cl <- makeCluster(10)
registerDoParallel(cl)
n_i <- c(80, 20)
# gamma <- c(1e-6, 0.01, 0.05, seq(0.1, 0.5, length.out = 5))
gamma <- 0
n_sims <- 1e4
results <- list()
idx <- 1
today <- as.numeric(format(Sys.time(), "%H%d%m%y"))
seed_start <- 10 * today
tic()
for(ii in 1:length(n_i)){
  for(jj in 1:length(gamma)){
    # Settings are (n, T, sigma, lambda_1, lambda_2, gamma)
    settings <- c(n_i[ii], 10, 1, 1, 1, gamma[jj])
    res_mat <- foreach(kk = 1:n_sims, .combine = rbind,
                       .errorhandling = "remove",
                       .packages = c("lmmstest", "lme4")) %dorng%{
                         c(do_one_sim(seed_start + kk, settings),
                           settings, seed_start + kk)
                       }
    results[[idx]] <- res_mat
    seed_start <- seed_start + n_sims
    idx <- idx + 1
    cat("Completed iteration", idx - 1, "of", length(n_i) * length(gamma) , "\n")
  }
}
toc()
stopCluster(cl)
out_mat <- do.call(rbind, results)
colnames(out_mat)[25:31] <- c("n", "r", "sig", "lam1_base", "lam2_base", "lam_scale", "seed")
my_path <- "~/"
saveRDS(out_mat, paste0(my_path, "lmm_sims_", today, ".Rds"))
