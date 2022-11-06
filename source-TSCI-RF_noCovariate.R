library(Matrix)
library(ranger)

# load('02_Bias correction/bias_correction_noCovariates_n5000_CF.RData')

# i <- 1
# Z <- cache.weight[i, ]$Z
# delta <- as.matrix(cache.weight[i, ]$Error[, 1]) 
# eps <- as.matrix(cache.weight[i, ]$Error[, 2]) 
# weight <- cache.weight[i, ]$weight
# # D
# f_1 <- function(z) {
# z + z^2 + 0.125 * z^4 - 25 / 12
# }
# D <- f_1(Z) + delta
# # Y
# Y <- D * gamma.true + eps
# # Arguments
# intercept=FALSE
# vio.space=NULL
# layer=TRUE
# split.prop=2/3
# Omega = weight
# num.trees=NULL
# mtry=NULL
# max.depth=NULL
# min.node.size=NULL
# str.thol=10
# alpha=0.05



# ==================================================================
TSCI.RF <- function(Y,D,Z,X=NULL,intercept=TRUE,vio.space=NULL,layer=TRUE,split.prop=2/3, Omega = NULL,
                    num.trees=NULL,mtry=NULL,max.depth=NULL,min.node.size=NULL,str.thol=10,alpha=0.05) {
  Y <- as.matrix(Y) 
  D <- as.matrix(D)
  Z <- as.matrix(Z)
  
  # constants
  n <- NROW(Y); p <- NCOL(Z)
  # default value for hyper-parameters
  if (is.null(num.trees)) num.trees <- 200
  if (is.null(mtry)) mtry <- seq(round(p/3),round(2*p/3),by=1)
  if (is.null(max.depth)) max.depth <- 0
  if (is.null(min.node.size)) min.node.size <- c(5,10,20)
  # define the vio.space if not specified
  if (is.null(vio.space)) {
    Q = 4
    vio.space <- matrix(NA,nrow(Z),0)
    for (q in 1:(Q-1)) {
      vio.space <- cbind(Z^q,vio.space)
    }
  } else {
    Q = NCOL(vio.space) + 1
  }
  # define the augmentation of covariates,
  # which is the combination of violation space and baseline covariates
  Cov.aug <- cbind(vio.space)
  n.A1 <- round(split.prop*n)
  A1.ind <- 1:n.A1
  
  if (is.null(Omega)) {
    # Treatment model fitting
    forest <- TSCI.RF.fit(D,Z,X,num.trees=num.trees,mtry=mtry,max.depth=max.depth,min.node.size=min.node.size,split.prop=split.prop)
    # Compute weight matrix
    weight <- TSCI.RF.weight(forest$nodes.A1)
  } else {
    weight = Omega
  }
  
  # Selection
  outputs <- TSCI.RF.Selection(Y,D,Cov.aug,A1.ind,weight=weight,Q=Q,intercept=intercept,layer=layer,str.thol=str.thol,alpha=alpha)
  return(outputs)
}


TSCI.RF.Selection <- function(Y, D, Cov.aug, A1.ind, weight, Q, intercept, layer, str.thol, alpha) {
  Y <- as.matrix(Y); D <- as.matrix(D); Cov.aug <- as.matrix(Cov.aug)
  # constants
  n <- length(Y); n.A1 <- length(A1.ind); r.aug <- NCOL(Cov.aug)
  Y.A1 <- Y[A1.ind]; D.A1 <- D[A1.ind]; Cov.aug.A1 <- Cov.aug[A1.ind,]
  # compute the representations
  Y.rep <- as.matrix(weight%*%Y.A1); D.rep <- as.matrix(weight%*%D.A1)
  Cov.rep <- as.matrix(weight%*%Cov.aug.A1)
  # the noise of treatment model
  delta.hat = D.A1 - D.rep
  SigmaSqD = mean(delta.hat^2)
  # save estimates for selection part
  names <- c(paste("RF-q",0:(Q-1),sep=""),paste("RF-Cor-q",0:(Q-1),sep=""))
  Coef.vec <- sd.vec <- rep(NA,2*Q)
  names(Coef.vec) <- names(sd.vec) <- names
  # IV strength test and signal strength test
  iv.str <- iv.thol <- rep(NA,Q)
  names(iv.str) <- names(iv.thol) <- paste("q",0:(Q-1),sep="")
  # the noise of outcome model
  eps.hat <- rep(list(NA),Q)
  # the numerator of variance
  explained.iv <- rep(NA,Q)
  names(explained.iv) <- paste("q",0:(Q-1),sep="")
  trace.T = explained.iv
  SigmaSqY = trace.T
  
  
  ### fixed violation space, compute necessary inputs of selection part
  # save D.resid for the computation of H and z.alpha
  D.resid <- RSS.vec <- rep(list(NA),Q)
  for (index in 1:Q) {
    q <- index-1
    if (q == Q-1) {
      if (intercept) {
        reg.rf <- lm(Y.rep~D.rep+Cov.rep)
        betaHat <- coef(reg.rf)[2]
      } else {
        reg.rf <- lm(Y.rep~D.rep+Cov.rep-1)
        betaHat <- coef(reg.rf)[1]
      }
      Coef.vec[index] <- betaHat
      eps.hat[[index]] = resid(lm(Y.A1-D.A1*betaHat~Cov.aug.A1))
      SigmaSqY[index] = mean(eps.hat[[index]]^2)
      stat.outputs <- TSCI.RF.stat(D.rep,Cov.rep,weight,n,eps.hat[[index]],delta.hat,str.thol=str.thol)
    } else if (q == 0) {
      if (intercept) {
        reg.rf <- lm(Y.rep~D.rep+0)
        betaHat <- coef(reg.rf)[2]
      } else {
        reg.rf <- lm(Y.rep~D.rep+0-1)
        betaHat <- coef(reg.rf)[1]
      }
      Coef.vec[index] <- betaHat
      eps.hat[[index]] = resid(lm(Y.A1-D.A1*betaHat~0))
      SigmaSqY[index] = mean(eps.hat[[index]]^2)
      stat.outputs <- TSCI.RF.stat(D.rep,Cov.rep[,-(1:(Q-1-q))],weight,n,eps.hat[[index]],delta.hat,str.thol=str.thol)
    } else {
      if (intercept) {
        reg.rf <- lm(Y.rep~D.rep+Cov.rep[,-(1:(Q-1-q))])
        betaHat <- coef(reg.rf)[2]
      } else {
        reg.rf <- lm(Y.rep~D.rep+Cov.rep[,-(1:(Q-1-q))]-1)
        betaHat <- coef(reg.rf)[1]
      }
      Coef.vec[index] <- betaHat
      eps.hat[[index]] = resid(lm(Y.A1-D.A1*betaHat~Cov.aug.A1[,-(1:(Q-1-q))]))
      SigmaSqY[index] = mean(eps.hat[[index]]^2)
      stat.outputs <- TSCI.RF.stat(D.rep,Cov.rep[,-(1:(Q-1-q))],weight,n,eps.hat[[index]],delta.hat,str.thol=str.thol)
    }
    # the statistics
    # Coef.vec[index+Q] <- stat.outputs$betaHat.cor
    sd.vec[index] <- stat.outputs$sd
    sd.vec[index+Q] <- stat.outputs$sd
    iv.str[index] <- stat.outputs$iv.str
    iv.thol[index] <- stat.outputs$iv.thol
    explained.iv[index] <- stat.outputs$explained.iv
    D.resid[[index]] <- stat.outputs$D.resid
    trace.T[index] = stat.outputs$trace.T
    RSS.vec[[index]] = stat.outputs$RSS.V
  }
  # Residual sum of squares of D.rep~Cov.rep
  D.RSS <- iv.str*SigmaSqD
  
  
  # # violation space selection
  # # all of the q below are from 0 to (Q-1), so use q+1 to index the columns
  # # Comparison and robust estimators
  # Coef.robust <- sd.robust <- rep(NA,4)
  # names(Coef.robust) <- names(sd.robust) <- c("RF-comp","RF-Cor-comp","RF-robust","RF-Cor-robust")
  # ivtest.vec <- (iv.str>=iv.thol)
  # run.OLS <- weak.iv <- FALSE
  # if (sum(ivtest.vec)==0) {
  #   warning("Weak IV, even if the IV is assumed to be valid; run OLS") # stop, output results of OLS
  #   run.OLS <- TRUE
  #   Qmax <- 1
  # } else {
  #   Qmax <- sum(ivtest.vec)-1
  #   if (Qmax==0) {
  #     warning("Weak IV, if the IV is invalid. We still test the IV invalidity.") # zijian: we need to rewrite this sentence...
  #     Qmax <- 1
  #     weak.iv = TRUE
  #   }
  # }
  # calculate Covariance at Qmax and compute bias correction
  # eps.Qmax = eps.hat[[Qmax+1]]
  # Coef.Qmax = rep(NA,Q)
  for (i in 1:Q) {
    # Coef.Qmax[i] = Coef.vec[i] - sum(RSS.vec[[i]]*delta.hat*eps.Qmax)/D.RSS[i]
    Coef.vec[i+Q] = Coef.vec[i] - sum(RSS.vec[[i]]*delta.hat*eps.hat[[i]])/D.RSS[i]
  }
  sd.vec[-(1:Q)] = sd.vec[1:Q]
  
  # ### Selection
  # # define comparison matrix
  # H <- beta.diff <- matrix(0,Qmax,Qmax)
  # # compute H matrix
  # for (q1 in 0:(Qmax-1)) {
  #   for (q2 in (q1+1):(Qmax)) {
  #     H[q1+1,q2] <- as.numeric(sum((weight%*%D.resid[[q1+1]])^2*eps.Qmax^2)/(D.RSS[q1+1]^2) + 
  #                                sum((weight%*%D.resid[[q2+1]])^2*eps.Qmax^2)/(D.RSS[q2+1]^2) - 
  #                                2*sum(eps.Qmax^2*(weight%*%D.resid[[q1+1]])*(weight%*%D.resid[[q2+1]]))/(D.RSS[q1+1]*D.RSS[q2+1])
  #     )
      
  #   }
  # }
  # # compute beta difference matrix, use Qmax
  # for (q in 0:(Qmax-1)) {
  #   beta.diff[q+1,(q+1):(Qmax)] <- abs(Coef.Qmax[q+1]-Coef.Qmax[(q+2):(Qmax+1)]) # use bias-corrected estimator
  # }
  # # bootstrap for the quantile of the differences
  # max.val <- rep(NA,300)
  # eps.Qmax.cent = eps.Qmax - mean(eps.Qmax)
  # for (i in 1:300) {
  #   diff.mat <- matrix(0,Qmax,Qmax)
  #   eps <- rep(NA,n.A1)
  #   for (j in 1:n.A1) {
  #     U.j = rnorm(1)
  #     eps[j] = eps.Qmax.cent[j]*U.j
  #   }
  #   eps.rep <- weight%*%eps
  #   for (q1 in 0:(Qmax-1)) {
  #     for (q2 in (q1+1):(Qmax)) {
  #       diff.mat[q1+1, q2] <- sum(D.resid[[q2+1]]*eps.rep)/(D.RSS[q2+1])-sum(D.resid[[q1+1]]*eps.rep)/(D.RSS[q1+1])
  #     }
  #   }
  #   diff.mat <- abs(diff.mat)/sqrt(H)
  #   max.val[i] <- max(diff.mat,na.rm = TRUE)
  # }
  # z.alpha <- 1.01*quantile(max.val,0.975)
  # diff.thol <- z.alpha*sqrt(H)
  # # comparison matrix
  # C.alpha <- ifelse(beta.diff<=diff.thol,0,1)
  
  # # layer selection or not
  # if (layer==TRUE) {
  #   # a vector indicating the selection of each layer
  #   sel.vec <- apply(C.alpha,1,sum)
  #   if (all(sel.vec != 0)) {
  #     q.comp = Qmax
  #   } else {
  #     q.comp = min(which(sel.vec==0))-1
  #   }
  # } else {
  #   ### What if Q = 3(q2) and Qmax = 1?
  #   sel.val <- C.alpha[1,Qmax]
  #   if (sel.val==1) {
  #     q.comp = Qmax
  #   } else {
  #     q.comp = 0
  #   }
  # }
  
  # ### invalidity of TSLS
  # if (q.comp>=1) {
  #   invalidity <- 1
  # } else {
  #   invalidity <- 0
  # }
  # q.robust <- min(q.comp+1, Qmax)
  # Coef.robust[1] <- Coef.vec[q.comp+1]
  # Coef.robust[2] <- Coef.vec[q.comp+Q+1]
  # Coef.robust[3] <- Coef.vec[q.robust+1]
  # Coef.robust[4] <- Coef.vec[q.robust+Q+1]
  # sd.robust[1] <- sd.vec[q.comp+1]
  # sd.robust[2] <- sd.vec[q.comp+Q+1]
  # sd.robust[3] <- sd.vec[q.robust+1]
  # sd.robust[4] <- sd.vec[q.robust+Q+1]
  # CI.robust = rbind(Coef.robust + qnorm(alpha/2)*sd.robust,Coef.robust + qnorm(1-alpha/2)*sd.robust)
  # rownames(CI.robust) = c("lower","upper")
  
  returnList = list(Coef.vec = Coef.vec,
                    sd.vec = sd.vec
                    # Coef.robust = Coef.robust,
                    # sd.robust = sd.robust,
                    # CI.robust = CI.robust,
                    # iv.str = iv.str, iv.thol = iv.thol,
                    # SigmaSqD = SigmaSqD,
                    # SigmaSqY = SigmaSqY,
                    # SigmaSqY.Qmax = mean(eps.Qmax^2),
                    # trace.T = trace.T,
                    # explained.iv = explained.iv,
                    # Qmax = Qmax,
                    # q.comp =q.comp, q.robust = q.robust,
                    # invalidity = invalidity,
                    # run.OLS = run.OLS,
                    # weak.iv = weak.iv
                    )
  returnList
}


TSCI.RF.stat <- function(D.rep, Cov.rep, weight, n, eps.hat, delta.hat, str.thol) {
  n.A1 <- length(D.rep); r.aug <- NCOL(Cov.rep)
  # compute the trace of T(V)
  # the trace of T matrix can be computed as RSS of each column of Omega on Cov.rep
  SigmaSqD = mean(delta.hat^2)
  RSS.V = rep(NA,n.A1)
  if (r.aug == 0) {
    for (j in 1:n.A1) {
      RSS.V[j] <- sum(resid(lm(weight[,j]~0))^2)
    }
    D.resid <- resid(lm(D.rep~0))
  } else {
    for (j in 1:n.A1) {
      RSS.V[j] <- sum(resid(lm(weight[,j]~Cov.rep))^2)
    }
    D.resid <- resid(lm(D.rep~Cov.rep))
  }
  
  trace.T = sum(RSS.V)
  D.rep2 <- as.matrix(weight%*%D.rep)
  # D.resid <- resid(lm(D.rep~Cov.rep))
  D.RSS <- sum(D.resid^2)
  iv.str <- D.RSS/SigmaSqD
  # this is the numerator of the variance of betaHat
  explained.iv <- as.numeric(t(D.resid)%*%weight%*%weight%*%D.resid)
  sd <- as.numeric(sqrt(sum(eps.hat^2*(weight%*%D.resid)^2))/D.RSS)
  
  
  
  ### standard errors of bias-corrected estimator
  # betaHat.cor <- betaHat - SigmaYD*trace.T/D.RSS
  # sd.cor <- sqrt((SigmaSqY*explained.iv+(SigmaSqD*SigmaSqY+SigmaYD^2)*(trace.T^2)/(n.A1-r.aug-1))/(D.RSS^2))
  
  # bootstrap for the threshold of IV strength test
  boot.vec <- rep(NA,300)
  delta.cent = delta.hat - mean(delta.hat)
  for (i in 1:300) {
    delta = rep(NA,n.A1)
    for (j in 1:n.A1) {
      U.j = rnorm(1)
      delta[j] = delta.cent[j]*U.j
    }
    
    delta.rep <- as.matrix(weight%*%delta)
    if (r.aug == 0) {
      delta.resid <- resid(lm(as.matrix(delta.rep)~0))
    } else {
      delta.resid <- resid(lm(as.matrix(delta.rep)~Cov.rep))
    }
    boot.vec[i] <- sum(delta.resid^2) + 2*sum(D.rep2*delta.resid)
  }
  iv.thol <- quantile(boot.vec,0.975)/SigmaSqD + max(2*trace.T, str.thol)
  # scale <- 1
  returnList <- list(
    # betaHat = betaHat,
    sd = sd,
    # betaHat.cor = betaHat.cor,
    # sd.cor = scale*sd.cor,
    D.resid = D.resid,
    iv.str = iv.str,
    iv.thol = iv.thol,
    explained.iv = explained.iv,
    trace.T = trace.T,
    RSS.V = RSS.V)
  returnList
}
