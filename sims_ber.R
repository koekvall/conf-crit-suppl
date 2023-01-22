#devtools::install_github("koekvall/lmmstest@paper-version")

# Load likelihood and auxiliary functions
source("bernoulli_funs.R")

do_one_sim <- function(seed, settings){

  set.seed(seed)

  n <- settings[1]
  r <- settings[2]
  mu <- settings[3]
  lam <- settings[4]
  num_nodes <- settings[5]
  dist <- ifelse(settings[6] == 1, "gamma", "normal")

  out_vec <- rep(0, 6)

  Y <- gen_mvber(n = n, mu = mu, lam = lam, r = r, dist = dist)

  # Proposed procedure
  tictoc::tic()
  s <- colSums(score(Y = Y, mu = mu, lam = lam, num_nodes = num_nodes,
                     dist = dist, modified = TRUE))
  I <- n * finf(mu = mu, lam = lam, r = r, num_nodes = num_nodes,
                dist = dist, modified = TRUE)
  our_stat <- crossprod(s, solve(I, s))
  timer <- tictoc::toc(quiet = TRUE)
  out_vec[1:2] <- c(our_stat, timer$toc - timer$tic)

  # Likelihood ratio procedure
  obj <- function(theta){
    -sum(log(f_marg(Y = Y, mu = theta[1], lam = theta[2], num_nodes = num_nodes, dist = dist)))
  }
  tictoc::tic()
  fit <- optim(c(mu, lam), fn = obj)
  lam_mle <- fit$par[2]
  mu_mle <- fit$par[1]
  lrt_stat <- 2 * sum(log(f_marg(Y = Y, mu = mu_mle, lam = lam_mle, num_nodes = num_nodes, dist = dist)/
                        f_marg(Y = Y, mu = mu, lam = lam, num_nodes = num_nodes, dist = dist)))
  timer <- tictoc::toc(quiet = TRUE)
  out_vec[3:4] <- c(lrt_stat, timer$toc - timer$tic)

  # Wald procedure (time is for evaluating test stat given MLEs)
  tictoc::tic()
  e <- c(mu_mle, lam_mle) - c(mu, lam)
  I_hat <- n * finf(mu = mu_mle, lam = lam_mle, r = r, num_nodes = num_nodes, dist = dist,
                    modified = FALSE)
  wald_stat <-  c(crossprod(e, I_hat %*% e))
  timer <- tictoc::toc(quiet = TRUE)
  out_vec[5:6] <- c(wald_stat, timer$toc - timer$tic)
  out_vec
}

# Do simulation
library(doParallel)
library(doRNG)
today <- as.numeric(format(Sys.time(), "%H%d%m%y"))
cl <- makeCluster(10)
registerDoParallel(cl)

n_sims <- 1e4
lam_vec <- c(0)

#dist <- "normal"

########################################
# To run with exponential random effects
dist <- "gamma"
#lam_vec <- sort(c(-lam_vec, lam_vec))
########################################


idx <- 1
seed_start <- 10 * today
out_list <- list()
for(ii in 1:length(lam_vec)){
  settings <- c(n = 80, r = 5, mu = 0.5, lam = lam_vec[ii], num_nodes = 10,
                dist = as.numeric(dist == "gamma")) # 1 means gamma, 0 means normal
  out_list[[ii]] <- foreach(kk = 1:n_sims,
                 .combine = rbind,
                 .packages = c("lme4", "tictoc"),
                 .errorhandling = "remove") %dorng%{
                   c(do_one_sim(seed_start + kk, settings), settings, seed_start + kk)
                 }

  cat("Completed simulation ", idx, "of ", 2 * length(lam_vec), "\n")
  idx <- idx + 1
  seed_start <- seed_start + n_sims
}
for(ii in 1:length(lam_vec)){
  settings <- c(n = 20, r = 5, mu = 0.5, lam = lam_vec[ii], num_nodes = 10,
                dist = as.numeric(dist == "gamma"))
  out_list[[ii + length(lam_vec)]] <- foreach(kk = 1:n_sims,
                            .combine = rbind,
                            .packages = c("lme4", "tictoc"),
                            .errorhandling = "remove") %dorng%{
                              c(do_one_sim(seed_start + kk, settings),
                                settings, seed_start + kk)
                            }
  cat("Completed simulation ", idx, "of ", 2 * length(lam_vec), "\n")
  idx <- idx + 1
  seed_start <- seed_start + n_sims
}
stopCluster(cl)
res_mat <- do.call(rbind, out_list)
colnames(res_mat)[1:6] <- c("stat_our", "time_our", "stat_lrt", "time_lrt",
                            "stat_wald", "time_wald")
colnames(res_mat)[12:13] <- c("dist", "seed")

my_path <- "~/"
saveRDS(res_mat, paste0(my_path, "ber_sims_", dist, today, ".Rds"))
