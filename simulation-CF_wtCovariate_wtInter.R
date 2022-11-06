library(parallel)
library(doParallel)
library(xtable)
source('source-weighting matrix.R')

# simulation setting
nsim <- 300
n <- 1000
p <- 5
inter <- TRUE
gamma.true <- c(0.1, 0.1, 0.2, 0.3, 0.4, 0.5)
eta.true <- rep(0.2, p)
rho1 <- 0.5
A1gen <- function(rho, p) {
  A1 <- matrix(0, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      A1[i, j] <- rho^(abs(i - j))
    }
  }
  A1
}


# Stage I: random forest (generate weighting matrix) -----------
registerDoParallel(12)
cache.weight <- foreach (i=1:nsim, .combine=rbind) %dopar% {
  library(Matrix)
  library(ranger)
  library(MASS)

    # X
  mu <- rep(0, p + 1)
  Cov <- (A1gen(rho1, p + 1))
  W.original <- mvrnorm(n, mu, Cov)
  W <- pnorm(W.original)
  X <- W[, -1]

    # Z
  Z <- W[, 1]
  Z <- 4*(Z-0.5)

    # Error (delta, epsilon)
  mu.error <- rep(0, 2)
  Cov.error <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  Error <- mvrnorm(n, mu.error, Cov.error)

    # D
  f_1 <- function(z) {
    z + z^2 + 0.125 * z^4 - 25 / 12
  }
  D <- f_1(Z) + Z*apply(X[, 1:3], 1, sum) - 0.3*apply(X, 1, sum)  + Error[, 1]

    # Y
  W <- X
  if (inter) {
    H <- apply(cbind(1, W), 2, function(x) x*D)
  } else {
    H <- as.matrix(D)
  }
  mu.X <- function(W, eta.true) {
    returntau = as.vector(W %*% eta.true)
    returntau
  }
  Y <- H %*% gamma.true + mu.X(W, eta.true) + Error[, 2]

  # the former part of function 'TSCI.RF.CF'
  Y <- as.matrix(Y)
  D <- as.matrix(D)
  Z <- as.matrix(Z)
  X <- as.matrix(X)
  n <- nrow(X)
  P <- ncol(X) + ncol(Z)

  # default value for hyper-parameters tuning
  num.trees <- 200
  mtry <- seq(round(P/3), round(2*P/3), by=1)
  max.depth <- 0
  min.node.size <- c(5, 10, 20)

  # Treatment model fitting
  forest <- TSCI.RF.fit(D, Z, X, num.trees = num.trees, mtry = mtry, max.depth = max.depth, min.node.size = min.node.size, split.prop = 2/3)
  # Compute weight matrix
  weight <- TSCI.RF.weight(forest$nodes.A1)
  # A1.ind <- forest$A1.ind
  # delta.hat <- D[A1.ind,] - weight%*%D[A1.ind,]

  cache <- list(X, Z, Error, weight)
  gc()
  return(cache)
}
colnames(cache.weight) <- c('X', 'Z', 'Error', 'weight')
gc()
stopImplicitCluster()

filename <- paste('weight-CF_wtCovariate_wtInter_n', n, '.RData', sep='')
save.image(filename)





# Stage II: control function part (estimate gamma) -----------------
filename <- paste('weight-CF_wtCovariate_wtInter_n', n, '.RData', sep='')
load(filename)


registerDoParallel(4)
stat <- foreach (i=1:nsim, .combine=rbind) %dopar% {
  library(Matrix)
  library(MASS)
  library(lmtest)

  X <- cache.weight[i, ]$X
  Z <- cache.weight[i, ]$Z
  delta <- as.matrix(cache.weight[i, ]$Error[, 1]) 
  eps <- as.matrix(cache.weight[i, ]$Error[, 2]) 
  weight <- as.matrix(cache.weight[i, ]$weight)

    # D
  f_1 <- function(z) {
    z + z^2 + 0.125 * z^4 - 25 / 12
  }
  D <- f_1(Z) + Z*apply(X[, 1:3], 1, sum) - 0.3*apply(X, 1, sum) + delta

    # Y
  W <- X
  if (inter) {
    H <- apply(cbind(1, W), 2, function(x) x*D)
  } else {
    H <- as.matrix(D)
  }
  mu.X <- function(W, eta.true) {
    returntau = as.vector(W %*% eta.true)
    returntau
  }
  Y <- H %*% gamma.true + mu.X(W, eta.true) + eps

  w.star <- c(1, rep(0.8, 5))
  true.val <- w.star %*% gamma.true
  split.prop <- 2/3
  n <- nrow(Y)
  n1 <- round(split.prop * n)
  A1.ind <- 1:n1
  D.A1 <- as.matrix(D[A1.ind, ])
  H.A1 <- as.matrix(H[A1.ind, ])
  W.A1 <- as.matrix(W[A1.ind, ])
  Y.A1 <- as.matrix(Y[A1.ind, ])
  delta.hat <- as.matrix(D.A1 - weight %*% D.A1)
  Delta <- diag(n1) - weight
  P <- ncol(H.A1)

  cf.result <- lm(Y.A1 ~ H.A1 + W.A1 + delta.hat - 1)
  eps.hat <- Y.A1 - H.A1 %*% coef(cf.result)[1:P] - W.A1 %*% coef(cf.result)[(P+1):(P+p)]

  # common statistics
  gamma.hat <- as.numeric(coef(cf.result)[1:P])
  bias <- as.numeric(gamma.hat - gamma.true)
  H.pro <- resid(lm(H.A1 ~ W.A1 + delta.hat))
  se <- sqrt(diag(solve(crossprod(H.pro)) * sum(eps.hat^2)/cf.result$df.residual))[1:P]
  c <- qt(0.975, cf.result$df.residual)
  CI <- c(gamma.hat - c*se, gamma.hat + c*se)
  CI.cover <- as.numeric(CI[1:P] <= gamma.true & gamma.true <= CI[(P+1):(2*P)])
  
  # values for estimating bias
  H.pro <- as.matrix(resid(lm(H.A1 ~ delta.hat + W.A1 - 1)))
  h <- solve(crossprod(H.pro))
  S <- list()
  for (j in 1:ncol(H.A1)) {
    if (j == 1) {
      S[[j]] <- diag(round(split.prop * n))
    } else {
      S[[j]] <- diag(W.A1[, j-1])
    }
  }
  W.A1.pro <- diag(n1) - W.A1 %*% solve(crossprod(W.A1)) %*% t(W.A1)
  temp <- W.A1.pro %*% Delta %*% D.A1
  deno <- crossprod(temp)
  rm(X, Z, delta, eps, D, H, Y, Y.A1, H.pro, temp)
  gc()

  # Eq29 hat
  B.29.hat <- rep(NA, ncol(H.A1))
  for (j in 1:ncol(H.A1)) {
    term1 <- t(D.A1) %*% S[[j]] %*% W.A1.pro %*% Delta %*% (D.A1) - t(delta.hat) %*% S[[j]] %*% W.A1.pro %*% Delta %*% (delta.hat)
    term2 <- t(delta.hat) %*% Delta %*% W.A1.pro %*% W.A1.pro %*% eps.hat
    B.29.hat[j] <- as.numeric(- term1 * term2 / deno)
  }
  bias.29.hat <- as.numeric(h %*% B.29.hat)

  # coverage
  CI.adj29Hat <- CI - as.numeric(bias.29.hat)
  CI.adj29Hat.cover <- as.numeric(CI.adj29Hat[1:P] <= gamma.true & gamma.true <= CI.adj29Hat[(P+1):(2*P)])
  se.val = sqrt(t(w.star) %*% vcov(cf.result)[1:6, 1:6] %*% w.star)
  val.hat = w.star %*% gamma.hat
  CI.val = c(val.hat - c*se.val, val.hat + c*se.val)
  CI.val.cover = as.numeric(CI.val[1] <= true.val & true.val <= CI.val[2])
  val.hat.adj29Hat = w.star %*% (gamma.hat - bias.29.hat)
  CI.val.adj29Hat = c(val.hat.adj29Hat - c*se.val, val.hat.adj29Hat + c*se.val)
  CI.val.adj29Hat.cover = as.numeric(CI.val.adj29Hat[1] <= true.val & true.val <= CI.val.adj29Hat[2])

  returnList <- c(gamma.hat = gamma.hat
                    , se = se
                    , bias = bias
                    , bias.29.hat = bias.29.hat
                    , CI.cover = CI.cover
                    , CI.adj29Hat.cover = CI.adj29Hat.cover
                    , CI.val.cover = CI.val.cover
                    , CI.val.adj29Hat.cover = CI.val.adj29Hat.cover
                  )
  gc()
  returnList
}
gc()
stopImplicitCluster()

filename <- paste0('simulation-CF_wtCovariate_wtInter_n', n, '.rds')
saveRDS(stat, filename)






# Analysis --------------------------------------------
filepath <- paste0('simulation-CF_wtCovariate_wtInter_n'
                  , c(1000, 3000, 5000)
                  , '.rds')
P = 6
bias.mat = se.mat = matrix(NA, 6, P)
cover.mat = matrix(NA, 6, P+1)
rownames(bias.mat) <- rownames(se.mat) <- rownames(cover.mat) <- paste0(rep(paste0('n=', c(1000, 3000, 5000)), each=2), rep(c('emp', 'est'), 2))
colnames(bias.mat) = colnames(se.mat) = paste0('gamma-', 1:P)
colnames(cover.mat) = c(paste0('gamma-', 1:P), 'w*gamma')
# colnames(cover.mat) = rep(c('Original coverage', 'Corrected coverage'), each=4)
for (i in 0:(length(filepath)-1)) {
  result <- readRDS(filepath[i+1])
  bias.mat[2*i+1, ] = colMeans(result)[13:(13+P-1)]
  bias.mat[2*i+2, ] = colMeans(result)[19:(19+P-1)]
  se.mat[2*i+1, ] = colMeans(result)[7:(7+P-1)]
  se.mat[2*i+2, ] = sqrt(colMeans((result[, 1:(1+P-1)] - matrix(rep(colMeans(result[, 1:(1+P-1)]), nrow(result)), nrow=nrow(result), byrow=TRUE))^2))
  cover.mat[2*i+1, 1:P] = colMeans(result)[25:(25+P-1)]
  cover.mat[2*i+1, -1:-P] = colMeans(result)[37]
  cover.mat[2*i+2, 1:P] = colMeans(result)[31:(31+P-1)]
  cover.mat[2*i+2, -1:-P] = colMeans(result)[38]
}
t(bias.mat)
t(se.mat)
t(cover.mat)

xtable(t(bias.mat), digits=4)
xtable(t(se.mat), digits=4)
xtable(t(cover.mat), digits=4)





