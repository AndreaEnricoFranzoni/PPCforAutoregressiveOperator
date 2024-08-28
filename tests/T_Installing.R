#to uploda the packages
#change here the directory
#setwd("/Users/andreafranzoni/Documents/Politecnico/Magistrale/PACS/PPCforAutoregressiveOperator")
#then 
#Rcpp::compileAttributes(".") 
# and then push the changes to modify the exports



library(Rcpp)
library(RcppEigen)



library(devtools)
devtools::install_github("AndreaEnricoFranzoni/PPCforAutoregressiveOperator", force = TRUE)

#function to evaluate mse on the predictions
mse = function(x1,x2)
{
  n = length(x1)
  err = numeric(n)
  
  for (i in 1:n) {
    err[i] = (x1[i] - x2[i])^2
  }
  
  return((sum(err)/n))
}


#PPC_KO version with R function
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
  
  ##add this
  #print(err.vec)
  
  alpha.opt <- alpha.vec[which.min(err.vec)]
  
  out <- list("alpha.opt" = alpha.opt,
              "err.vec"   = err.vec,
              "alpha.vec" = alpha.vec)
  return(out)
}
