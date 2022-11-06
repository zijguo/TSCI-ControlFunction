
library(ranger)

# fit random forest
TSCI.RF.fit <- function(D, Z, X, num.trees, mtry, max.depth, min.node.size, split.prop, MSE.thol=1e12, forest.save = TRUE) {
  W <- as.matrix(cbind(Z, X))
  D <- as.matrix(D)
  n <- NROW(W)
  p <- NCOL(W)
  Data <- data.frame(cbind(D, W))
  names(Data) <- c("D", paste("W", 1:p, sep = ""))

  # grid search
  params.grid <- expand.grid(
    num.trees = num.trees,
    mtry = mtry,
    max.depth = max.depth,
    min.node.size = min.node.size
  )

  # split the data into two parts A1 and A2
  n.A1 <- round(split.prop * n)
  A1.ind <- 1:n.A1
  Data.A1 <- Data[A1.ind, ]
  Data.A2 <- Data[-A1.ind, ]

  # use A2 to build the random forest
  forest.A2 <- NULL
  MSE.oob.A2 <- MSE.thol
  params.A2 <- NULL
  ### use oob error to do hyper-parameter tuning
  for (i in 1:nrow(params.grid)) {
    temp.A2 <- ranger(D ~ .,
      data = Data.A2,
      num.trees = params.grid$num.trees[i],
      mtry = params.grid$mtry[i],
      max.depth = params.grid$max.depth[i],
      min.node.size = params.grid$min.node.size[i],
      importance = "impurity"
    )
    if (temp.A2$prediction.error <= MSE.oob.A2) {
      forest.A2 <- temp.A2
      params.A2 <- params.grid[i, ]
      MSE.oob.A2 <- temp.A2$prediction.error
    }
  }

  # leaf nodes information of A1 on the Random Forest built on A2
  nodes.A1 <- predict(forest.A2, data=Data.A1, type="terminalNodes")$predictions
  returnList <- list(
    forest.A2 = forest.A2,
    params.A2 = params.A2,
    A1.ind = A1.ind,
    nodes.A1 = nodes.A1,
    MSE.oob.A2 = MSE.oob.A2
  )
  if (!forest.save) returnList <- returnList[-1]
  returnList
}

# generate weight matrix
TSCI.RF.weight <- function(nodes) {
  n.A1 <- NROW(nodes)
  num.trees <- NCOL(nodes)
  out.weight <- matrix(0, n.A1, n.A1)
  for (j in 1:num.trees) {
    weight.mat <- matrix(0, n.A1, n.A1) # weight matrix for single tree
    unique.nodes <- unique(nodes[, j])
    for (i in 1:length(unique.nodes)) {
      ind <- nodes[, j] == unique.nodes[i] # indices of samples in the node
      num.samples <- sum(ind) # number of samples in the node
      w <- 1 / (num.samples - 1) # weight, to remove self-prediction
      weight.vec <- ifelse(ind, yes = w, no = 0)
      weight.mat[ind, ] <- matrix(rep(weight.vec, num.samples), num.samples, byrow = T) / num.trees
    }
    diag(weight.mat) <- 0 # remove self prediction
    out.weight <- out.weight + weight.mat
  }
  out.weight <- Matrix(out.weight, sparse = T) # sparse matrix to save memory
  out.weight
}
