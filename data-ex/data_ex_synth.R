source("bernoulli_funs.R")
set.seed(33)

# Settings
n <- 80
r <- 5
num_nodes <- 10
mu0 <- 0.5
lam0 <- -0.5

# Generate data
Y <- gen_mvber(n = n, mu = mu0, lam = lam0, r = r, dist = "gamma")

# Proposed test-statistic
compute_stat <- function(theta, idx = 1:length(theta)){
   s <- colSums(score(Y, theta[1], theta[2], num_nodes = num_nodes,
                      dist = "gamma",
                      modified = TRUE))[idx]
   I <- finf(theta[1], theta[2], r = r, num_nodes = num_nodes,
             dist = "gamma", modified = TRUE)
   if(length(idx) != 2){
     out <- (s^2 / (I[idx, idx] - I[idx, -idx]^2 / I[-idx, -idx])) / n
   } else{
     out <- c(crossprod(s, solve(I, s))) / n
   }
   out
}

# Maximize crude approximation to get idea of where to center region
neg_ll_approx <- function(theta){
  -sum(log(f_marg(Y = Y, mu = theta[1], lam = theta[2], num_nodes = 1,
                  dist = "gamma")))
}
tictoc::tic()
optim(c(0, 0), neg_ll_approx)
tictoc::toc()

# Check componentwise statistics to get an idea of how large to make region
# For fixed effect
compute_stat(c(-1, -0.044), idx = 1)
compute_stat(c(1, -0.044), idx = 1)
compute_stat(c(0, -0.044), idx = 1)
# For scale parameter
compute_stat(c(0.55, -1), idx = 2)
compute_stat(c(0.55, 1), idx = 2)
compute_stat(c(0.55, 0), idx = 2)

# Construct confidence region
mu_vals <- seq(0, 1, length.out = 50)
lam_vals <- seq(-1, 1, length.out = 50)
stat_mat <- matrix(0, length(mu_vals), length(lam_vals))
covered_theta <- list()
idx <- 1
tictoc::tic()
for(ii in 1:length(mu_vals)){
 cat(ii, " of ", length(mu_vals), "\n")
 for(jj in 1:length(lam_vals)){
   stat_mat[ii, jj] <- compute_stat(c(mu_vals[ii], lam_vals[jj]))
   if(stat_mat[ii, jj] <= qchisq(0.95, 2)){
     covered_theta[[idx]] <- c(mu_vals[ii], lam_vals[jj])
     idx <- idx + 1
   }
 }
}
tictoc::toc()
colnames(stat_mat) <- round(lam_vals, 2)
rownames(stat_mat) <- round(mu_vals, 2)

# Plot confidence region
out_dir <- "~/" # end in "/"
PDF <- TRUE
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2", "#D55E00", "#CC79A7")
if(PDF) pdf(paste(out_dir, "fig_data_ex_synth.pdf", sep = ""), width = 9, height = 3)
par(cex.axis = 1.1, cex.lab = 1.1, cex = 1.1)
par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.5,1.1))

contour(x = mu_vals,
       y = lam_vals,
       z = stat_mat,
       levels = c(0, round(qchisq(c(0.8, 0.9, 0.95), df = 2), digits = 2)),
       drawlabels = TRUE,
       labcex = 1,
       labels = c(0, 0.8, 0.9, 0.95),
       xlab = bquote(~ psi),
       ylab = bquote(~ lambda))
dev.off()

# Image of confidence region
possible_mean <- rep(0, length(covered_theta))
for(ii in 1:length(covered_theta)){
  possible_mean[ii] <- f_marg(matrix(1, 1, 1),
                              covered_theta[[ii]][1],
                              covered_theta[[ii]][2],
                              num_nodes = num_nodes,
                              dist = "gamma")
}
range(possible_mean)

# True success probability
f_marg(matrix(1, 1, 1), mu0, lam0, num_nodes = num_nodes, dist = "gamma")
