#' ESTIMATION OF THE LAG-1 AUTOREGRESSIVE OPERATOR

## MAIN FUNCTIONS --------------------------------------------------------------
## Principal Predictive Components - Regularized strategy
PPCdiscretized_KO <- function(X, p.thresh=0.95, alpha=NULL, K=NULL)
{
  n <- dim(X)[2]
  TT <- dim(X)[1] # number of measurements points (of sampling points on the functional domain)
  # row.means <- rowMeans(X)
  # X <- X - row.means
  
  ## Estimation of the covariance and its inverse over the subspace spanned by the first
  ## K eigenvectors
  autocov <- tcrossprod(X)/n
  crosscov <- tcrossprod(X[,2:n], X[,1:(n-1)])/(n-1)
  
  if(is.null(alpha)){
    alpha <- 0.75 
  }
  
  Gamma1_2 <- t(crosscov) %*% crosscov
  Gamma0 <- autocov + alpha*diag(TT)
  
  eigendecB <- eigen(Gamma0)
  lambda.alpha <- eigendecB$values
  W.alpha <- eigendecB$vectors
  p <- head(which(cumsum(lambda.alpha)/sum(lambda.alpha)>=p.thresh),n=1)
  if(p>1){L.alpha <- diag(1/sqrt(lambda.alpha[1:p]))}else{L.alpha <- 1/sqrt(lambda.alpha)[1]}
  W.alpha <- matrix(W.alpha[,1:p], nrow=dim(W.alpha)[1], ncol=p)

  sqGamma0.inv <- W.alpha %*% L.alpha %*% t(W.alpha)

  Phi <- sqGamma0.inv %*% Gamma1_2 %*% sqGamma0.inv
  
  eigendec <- eigen(Phi)
  V1 <- eigendec$vectors
  D1 <- eigendec$values
  V1 <- V1[,order(D1, decreasing=TRUE)]
  D1 <- sort(D1, decreasing=TRUE)
  
  if(is.null(K)){
    K <- p 
  }
  
  Vhat1 <- V1[,1:K]
  Dhat1 <- D1[1:K]
  
  VVT <- Vhat1 %*% t(Vhat1)
  Psi.hat <- crosscov %*% sqGamma0.inv %*% VVT %*% sqGamma0.inv
  
  return(Psi.hat)
  
}



## Principal Predictive Components - GEP Strategy
PPCdiscretized <- function(X, alpha=NULL, K=NULL)
{
  n <- dim(X)[2]
  TT <- dim(X)[1] # number of measurements points (of sampling points on the functional domain)
  # row.means <- rowMeans(X)
  # X <- X - row.means
  
  ## Estimation of the covariance and its inverse over the subspace spanned by the first
  ## K eigenvectors
  autocov <- tcrossprod(X)/n
  crosscov <- tcrossprod(X[,2:n], X[,1:(n-1)])/(n-1)
  
  Gamma1_2 <- t(crosscov) %*% crosscov
  Gamma0 <- autocov
  
  if(is.null(alpha)){
    alpha.vec <- seq(1, 15, by=1)
    alpha.sel <- dalpha_select(alpha.vec, XX=X, K=TT)
    alpha.step1 <- alpha.sel$alpha.opt
    
    alpha.vec <- seq(max(0.1, alpha.step1-1), alpha.step1+1, by=0.1)
    alpha.sel <- dalpha_select(alpha.vec, XX=X, K=TT)
    alpha.step2 <- alpha.sel$alpha.opt
    
    alpha.vec <- seq(max(0.01,alpha.step2-0.1), alpha.step2+0.1, by=0.01)
    alpha.sel <- dalpha_select(alpha.vec, XX=X, K=TT)
    alpha <- alpha.sel$alpha.opt 
  }
  
  geigendec <- geigen::geigen(Gamma1_2, Gamma0 + alpha*diag(TT))
  V1 <- geigendec$vectors
  D1 <- geigendec$values
  V1 <- V1[,order(D1, decreasing=TRUE)]
  D1 <- sort(D1, decreasing=TRUE)
  
  if(is.null(K)){
    nppc.sel <- dK_select(XX=X, alpha)
    K <- nppc.sel$nppc.opt 
  }
  
  Vhat1 <- V1[,1:K]
  Dhat1 <- D1[1:K]
  
  ## Allora i principal predictive factors sono le colonne di W = sqrt(Gamma0) %*% Vhat1
  ## Invece non ho pensato ai loadings
  
  Psi.hat <- crosscov %*% Vhat1 %*% t(Vhat1)
  return(Psi.hat)
  
}

## Yule-Walker predictor: Psi = Gamma0 %*% Gamma^(-1)
YWdiscretized <- function(X)
{
  n <- dim(X)[2]
  TT <- dim(X)[1] # number of measurements points (of sampling points on the functional domain)
  # row.means <- rowMeans(X)
  # X <- X - row.means
  
  ## Estimation of the covariance and its inverse over the subspace spanned by the first
  ## nharm eigenvectors
  autocov <- tcrossprod(X)/n
  eigendec <- eigen(autocov)
  
  nharm <- dharm_select(XX=X)$nharm.opt
  
  D <- eigendec$values[1:nharm]
  Vhat <- eigendec$vectors[,1:nharm]
  D.inv <- 1/D
  # Gamma1_2 <- matrix(data=rep(1,TT), ncol=1) %*% D.inv
  # vlambda <- Vhat * Gamma1_2
  # acovinv <- vlambda %*% t(Vhat)
  if(nharm > 1){
    acovinv <- Vhat %*% diag(D.inv) %*% t(Vhat)
  } else {
    acovinv <- D.inv * tcrossprod(Vhat)
  }
  
  ## Estimation of the cross covariance and the autocorrelation
  crosscov <- tcrossprod(X[,2:n], X[,1:(n-1)])/(n-1)
  autocor <- tcrossprod(Vhat) %*% crosscov %*% acovinv
  
  return(autocor)
  
}

## Estimated Kernel as in Didericksen, Kokoszka and Zhang (2012)
EKdiscretized <- function(X, p.thresh=0.95, p = NULL)
{
  n <- dim(X)[2]
  TT <- dim(X)[1] # number of measurements points (of sampling points on the functional domain)
  # row.means <- rowMeans(X)
  # X <- X - row.means
  
  ## Estimation of the covariance and its inverse over the subspace spanned by the first
  ## K eigenvectors
  autocov <- tcrossprod(X)/n
  crosscov <- tcrossprod(X[,2:n], X[,1:(n-1)])/(n-1)
  
  alpha <- 0
  
  Gamma0 <- autocov + alpha*diag(TT)
  
  eigendecB <- eigen(Gamma0)
  lambda.alpha <- eigendecB$values
  W.alpha <- eigendecB$vectors
  
  if( is.null(p) )
  {
    p <- head(which(cumsum(lambda.alpha)/sum(lambda.alpha)>=p.thresh),n=1)
  }
  
  if(p>1){L.alpha <- diag(1/sqrt(lambda.alpha[1:p]))}else{1/sqrt(lambda.alpha)[1:p]}
  W.alpha <- matrix(W.alpha[,1:p], nrow=dim(W.alpha)[1], ncol=p)
  
  Psi.mat <- 1/(TT^2) * L.alpha %*% t(W.alpha) %*% crosscov %*% W.alpha
  Psi.hat <- 1/TT*W.alpha %*% Psi.mat %*% t(W.alpha)
  
  return(Psi.hat)
  
}

## Estimated Kernel Improved as in Didericksen, Kokoszka and Zhang (2012)
EKimproved <- function(X, p.thresh=0.95, p=NULL)
{
  n <- dim(X)[2]
  TT <- dim(X)[1] # number of measurements points (of sampling points on the functional domain)
  # row.means <- rowMeans(X)
  # X <- X - row.means
  
  ## Estimation of the covariance and its inverse over the subspace spanned by the first
  ## K eigenvectors
  autocov <- tcrossprod(X)/n
  crosscov <- tcrossprod(X[,2:n], X[,1:(n-1)])/(n-1)
  
  alpha <- 0
  
  Gamma0 <- autocov + alpha*diag(TT)
  
  eigendecB <- eigen(Gamma0)
  lambda.alpha <- eigendecB$values
  lambda.alpha[3:length(lambda.alpha)] <- lambda.alpha[3:length(lambda.alpha)] + 1.5*(lambda.alpha[1]+lambda.alpha[2])
  W.alpha <- eigendecB$vectors
  
  if( is.null(p) )
  {
    p <- head(which(cumsum(lambda.alpha)/sum(lambda.alpha)>=p.thresh),n=1)
  }
  
  if(p>1){L.alpha <- diag(1/sqrt(lambda.alpha[1:p]))}else{L.alpha <- 1/sqrt(lambda.alpha)[1:p]}
  W.alpha <- matrix(W.alpha[,1:p], nrow=dim(W.alpha)[1], ncol=p)
  
  Psi.mat <- 1/(TT^2) * L.alpha %*% t(W.alpha) %*% crosscov %*% W.alpha
  Psi.hat <- 1/TT*W.alpha %*% Psi.mat %*% t(W.alpha)
  
  return(Psi.hat)
  
}

## UTILITIES -------------------------------------------------------------------
#' Hyperparameter selection in an increasing-window fashion

dalpha_select <- function(alpha.vec, XX, K=1)
{
  n <- dim(XX)[2]
  TT <- dim(XX)[1]
  
  # Errors' vector
  err.vec <- numeric(length(alpha.vec))
  
  for(l in 1:length(alpha.vec)) # l=1
  {
    alpha <- alpha.vec[l]

    err.cv <- numeric(0)
    
    for(b in ceiling(n/2):(n-1))
    {
      X.train  <- XX[,1:b]
      n.train  <- dim(X.train)[2]
      autocov  <- tcrossprod(X.train)/n.train
      crosscov <- tcrossprod(X.train[,2:n.train], X.train[,1:(n.train-1)])/(n.train-1)
      
      ## Define and solve the generalized eigenvalues problem
      Gamma1_2 <- t(crosscov) %*% crosscov
      Gamma0 <- autocov
      
      geigendec <- geigen::geigen(Gamma1_2, Gamma0 + alpha*diag(TT))
      V1 <- geigendec$vectors
      D1 <- geigendec$values
      V1 <- V1[,order(D1, decreasing=TRUE)]
      D1 <- sort(D1, decreasing=TRUE)
      
      Vhat1 <- V1[,1:K]
      Dhat1 <- D1[1:K]
      
      Psi.hat <- crosscov %*% Vhat1 %*% t(Vhat1)
      
      ## Evaluate the error on the test set with alpha
      ## Psi
      y.true <- XX[,(b+1)]
      ## Psi hat
      y.hat <- Psi.hat %*% XX[,b]
      
      err <- y.true - y.hat
      err.cv <- c(err.cv,mean(err^2))
      
    }
    err.vec[l] <- mean(err.cv)
  }

  alpha.opt <- alpha.vec[which.min(err.vec)]
  
  out <- list("alpha.opt" = alpha.opt,
              "err.vec"   = err.vec,
              "alpha.vec" = alpha.vec)
  return(out)
}

dK_select <- function(XX, alpha)
{
  n <- dim(XX)[2]
  TT <- dim(XX)[1]
  n.train <- round(n/2)
  
  X.train <- XX[,1:n.train]
  autocov <- tcrossprod(X.train)/n.train
  crosscov <- tcrossprod(X.train[,2:n.train], X.train[,1:(n.train-1)])/(n.train-1)
  
  ## Define and solve the generalized eigenvalues problem
  Gamma1_2 <- t(crosscov) %*% crosscov
  Gamma0 <- autocov
  geigendec <- geigen::geigen(Gamma1_2, Gamma0 + alpha*diag(TT))
  V1 <- geigendec$vectors
  D1 <- geigendec$values
  V1 <- V1[,order(D1, decreasing=TRUE)]
  D1 <- sort(D1, decreasing=TRUE)
  
  nppc.vec <- which(D1>=1e-4)
  nppc.vec <- which(cumsum(D1)/sum(D1) < 1)
  err.vec <- numeric(length(nppc.vec))
  
  for(l in 1:length(nppc.vec)) # l=1
  {
    K <- nppc.vec[l]
    
    Vhat1 <- V1[,1:K]
    Dhat1 <- D1[1:K]
    
    Psi.hat <- crosscov %*% Vhat1 %*% t(Vhat1)
    
    ## Evaluate the error on the test set with alpha
    ## Psi
    y.true <- XX[,(n.train+2):n]
    ## Psi hat
    y.hat <- Psi.hat %*% XX[,(n.train+1):(n-1)]
    
    err <- y.true - y.hat
    err.vec[l] <- sqrt(sum(err^2))
    
  }
  
  nppc.opt <- nppc.vec[which.min(err.vec)]
  
  out <- list("nppc.opt" = nppc.opt,
              "err.vec"   = err.vec,
              "nppc.vec" = nppc.vec)
  return(out)
  
}

dharm_select <- function(XX)
{
  n <- dim(XX)[2]
  n.train <- round(n/2)
  X.train <- XX[,1:n.train]
  
  autocov <- tcrossprod(X.train)/n.train
  crosscov <- tcrossprod(X.train[,2:n.train], X.train[,1:(n.train-1)])/(n.train-1)
  eigendec <- eigen(autocov)
  
  nharm.vec <- which(eigendec$values >= 1e-4)
  err.vec <- numeric(length(nharm.vec))
  for(l in 1:length(nharm.vec))
  {
    nharm <- nharm.vec[l]
    
    D <- eigendec$values[1:nharm]
    Vhat <- eigendec$vectors[,1:nharm]
    D.inv <- 1/D
    
    if(nharm > 1){
      acovinv <- Vhat %*% diag(D.inv) %*% t(Vhat)
    } else {
      acovinv <- D.inv * tcrossprod(Vhat)
    }
    
    ## Estimation of the cross covariance and the autocorrelation
    Psi.hat <- tcrossprod(Vhat) %*% crosscov %*% acovinv
    
    ## Evaluate the error on the test set with nharm
    ## Psi
    y.true <- XX[,(n.train+2):n]
    ## Psi hat
    y.hat <- Psi.hat %*% XX[,(n.train+1):(n-1)]
    
    err <- y.true - y.hat
    err.vec[l] <- sqrt(sum(err^2))
    
  }
  
  nharm.opt <- nharm.vec[which.min(err.vec)]
  
  out <- list("nharm.opt" = nharm.opt,
              "err.vec"   = err.vec,
              "nharm.vec" = nharm.vec)
  return(out)
}
