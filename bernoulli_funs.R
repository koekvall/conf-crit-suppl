logist <- function(t){1 / (1 + exp(-t))}

f_cond <- function(Y, w, mu, lam){
  # Y and w should be be n x r and n x 1, respectively
  # mu and lam should be scalars
  r <- ncol(Y)
  y <- as.vector(Y)
  P <- matrix(dbinom(y,
                     size = 1,
                     prob = logist(mu + lam * rep(w, r))),
              ncol = r)
  
  # Return
  apply(P, 1, prod)
}

f_marg <- function(Y, mu, lam, num_nodes, dist){
  # Y should be n x r matrix; mu and lam scalar
  quad <- statmod::gauss.quad.prob(n = num_nodes, dist = dist)
  quad$nodes <- quad$nodes - ifelse(dist == "gamma", 1, 0)
  
  # Return
  apply(Y, 1, function(y){
    sum(quad$weights * f_cond(Y = matrix(rep(y, each = num_nodes),
                                         nrow = num_nodes),
                              w = quad$nodes,
                              mu = mu,
                              lam = lam))
  })
}

score <- function(Y, mu, lam, num_nodes = 10, dist, modified = FALSE,
                  tol = 1e-12){
  # Y should be n x r matrix
  # mu should be scalar
  n <- nrow(Y); r <- ncol(Y)
  quad <- statmod::gauss.quad.prob(n = num_nodes, dist = dist)
  quad$nodes <- quad$nodes - ifelse(dist == "gamma", 1, 0)
  score_mu <- apply(Y, 1, function(y){
    -sum(quad$weights * (r * logist(mu + quad$nodes * lam) - sum(y)) *
           f_cond(Y = matrix(rep(y, each = num_nodes),
                             nrow = num_nodes),
                  w = quad$nodes,
                  mu = mu,
                  lam = lam))
  })
  if((abs(lam) > tol) | (!modified)){
    score_lam <- apply(Y, 1, function(y){
      -sum(quad$weights * quad$nodes *
             (r * logist(mu + quad$nodes * lam) - sum(y)) *
             f_cond(Y = matrix(rep(y, each = num_nodes), nrow = num_nodes),
                    w = quad$nodes, mu = mu, lam = lam))
    })
  } else{
    score_lam <- apply(Y, 1, function(y){
      prob <- logist(mu + quad$nodes * lam)
      -sum(quad$weights * quad$nodes^2 *
             ((r * prob - sum(y))^2 - r * prob * (1 - prob)) *
             f_cond(Y = matrix(rep(y, each = num_nodes),
                               nrow = num_nodes), w = quad$nodes, mu = mu,
                    lam = lam))
    })
  }
  # Return
  cbind(score_mu, score_lam) * (1 / f_marg(Y, mu, lam, num_nodes, dist = dist))
}

finf <- function(mu, lam, r, num_nodes, dist, modified = FALSE, tol = 1e-12)
{
  Y <- matrix(rep(0, r), nrow = 1)
  I <- tcrossprod(c(score(Y = Y, mu = mu, lam = lam,
                          num_nodes = num_nodes, dist = dist,
                          modified = modified, tol = tol))) *
    f_marg(Y = Y, mu = mu, lam = lam, num_nodes = 10, dist = dist)
  for(ii in 1:r){
    Y[1, ii] <- 1
    I <- I + choose(r, ii) * tcrossprod(c(score(Y = Y, mu = mu, lam = lam,
                                                num_nodes = num_nodes,
                                                dist = dist,
                                                modified = modified,
                                                tol = tol))) *
      f_marg(Y = Y, mu = mu, lam = lam, num_nodes = 10, dist = dist)
  }
  
  #Return
  I
}

gen_mvber <- function(n, mu, lam, r, dist){
  if(dist == "normal"){
    w <- rnorm(n)
  } else{
    w <- rexp(n) - 1
  }
  lin_pred <- matrix(lam * rep(w, r), ncol = r) + mu
  # Return
  matrix(rbinom(n = r * n, size = 1, prob = logist(as.vector(lin_pred))), ncol = r)
}

