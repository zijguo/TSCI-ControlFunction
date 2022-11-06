library(parallel)
library(doParallel)
library(xtable)
source('source-weighting matrix.R')

# simulation setting
nsim <- 300
n <- 1000
p <- 5
gamma.true <- 0.2
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


# Stage I: random forest (generate weighting matrix) --------------------------
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
  D <- f_1(Z) + Error[, 1]
  
  # Y
  Y <- D * gamma.true + Error[, 2]
  
  # the former part of function 'TSCI.RF.CF'
  Y <- as.matrix(Y)
  D <- as.matrix(D)
  Z <- as.matrix(Z)
  n <- nrow(Y)
  
  # default value for hyper-parameters tuning
  num.trees <- 200
  mtry <- 1
  max.depth <- 0
  min.node.size <- c(5, 10, 20)
  
  # Treatment model fitting
  forest <- TSCI.RF.fit(D, Z, X=NULL, num.trees = num.trees, mtry = mtry, max.depth = max.depth, min.node.size = min.node.size, split.prop = 2/3)
  # Compute weight matrix
  weight <- TSCI.RF.weight(forest$nodes.A1)
  
  cache <- list(Z, Error, weight)
  gc()
  return(cache)
}
colnames(cache.weight) <- c('Z', 'Error', 'weight')
gc()
stopImplicitCluster()

filename <- paste('weight-CF_noCovariates_n', n, '.RData', sep='')
save.image(filename)







# Stage II: control function part (estimate gamma) -------------------
filename <- paste('weight-CF_noCovariates_n', n, '.RData', sep='')
load(filename)


registerDoParallel(12)
stat <- foreach (i=1:nsim, .combine=rbind) %dopar% {
  library(Matrix)
  library(MASS)
  library(lmtest)
  
  Z <- cache.weight[i, ]$Z
  delta <- as.matrix(cache.weight[i, ]$Error[, 1]) 
  eps <- as.matrix(cache.weight[i, ]$Error[, 2]) 
  weight <- as.matrix(cache.weight[i, ]$weight)
  
  # D
  f_1 <- function(z) {
    z + z^2 + 0.125 * z^4 - 25 / 12
  }
  D <- f_1(Z) + delta
  
  # Y
  Y <- D * gamma.true + eps
  
  split.prop <- 2/3
  n <- nrow(Y)
  n1 <- round(split.prop * n)
  A1.ind <- 1:n1
  Z.A1 <- as.matrix(Z[A1.ind, ])
  eps.A1 <- as.matrix(eps[A1.ind, ])
  D.A1 <- as.matrix(D[A1.ind, ])
  Y.A1 <- as.matrix(Y[A1.ind, ])
  f.hat <- weight %*% D.A1
  delta.hat <- D.A1 - f.hat
  Delta <- diag(length(A1.ind)) - weight

  pro.matrix.comp <- diag(length(A1.ind)) - (tcrossprod(delta.hat)) / as.numeric(crossprod(delta.hat))
  Y.pro <- as.matrix(pro.matrix.comp%*%Y.A1)
  D.pro <- as.matrix(pro.matrix.comp%*%D.A1)
  cf.result <- lm(Y.A1 ~ D.A1 + delta.hat - 1)
  eps.hat <- Y.A1 - D.A1*coef(cf.result)[1]
  
  # common statistics
  gamma.hat <- as.numeric(coef(cf.result)[1])
  bias <- gamma.hat - gamma.true
  se <- sqrt(solve(crossprod(D.pro)) * sum(eps.hat^2)/cf.result$df.residual)
  c <- qt(0.975, cf.result$df.residual)
  CI <- c(gamma.hat - c * se, gamma.hat + c * se)
  CI.cover <- (CI[1] <= gamma.true) & (gamma.true <= CI[2])

  returnList <- c(gamma.hat = gamma.hat
                    , bias = bias
                    , se = se
                    , CI.cover = CI.cover
  )
  gc()
  returnList
}
gc()
stopImplicitCluster()

filename <- paste0('simulation-CF_noCovariates_n', n, '.rds')
saveRDS(stat, filename)






# Analysis -----------------------------------------
files = paste0('simulation-CF_noCovariates_n'
                , c(1000, 3000, 5000)
                , '.rds')

bias.mat = cover.mat = matrix(NA, 3, 1)
se.mat = matrix(NA, 3, 2)
rownames(bias.mat) = rownames(se.mat) = rownames(cover.mat) = paste0('n=', c(1000, 3000, 5000))
colnames(bias.mat) = c('Empirical bias')
colnames(se.mat) = c('Calculated SE', 'Empirical SE')
colnames(cover.mat) = c('Empirical coverage')
for (i in 1:length(files)) {
  stat = readRDS(files[i])
  bias.mat[i, ] = colMeans(stat)[2]
  se.mat[i, 1] = colMeans(stat)[3]
  se.mat[i, 2] = sqrt(mean((stat[, 1] - colMeans(stat)[1])^2))
  cover.mat[i, ] = colMeans(stat)[4]
}
t(bias.mat)
t(se.mat)
t(cover.mat)

xtable(t(bias.mat), digits=4)
xtable(t(se.mat), digits=4)
xtable(t(cover.mat), digits=4)