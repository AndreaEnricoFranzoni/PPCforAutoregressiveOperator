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



#to evaluate mse of prediction
mse = function(x1,x2)
{
  n = length(x1)
  err = numeric(n)
  
  for (i in 1:n) {
    err[i] = (x1[i] - x2[i])^2
  }
  
  return((sum(err)/n))
}


set.seed(230300)
library(Rcpp)
library(RcppEigen)

#change here the directory
setwd("/Users/andreafranzoni/Documents/Politecnico/Magistrale/PACS/PPCforAutoregressiveOperator")

#installing the package
#Rcpp::compileAttributes("PPCKO") 
devtools::install("PPCKO")

#leave these parameters like that for now
cv = "NoCV"
threshold=0.95
alpha=0.75
n_disc=100



#simple test
data_test = matrix(c(1,4,7,10,6,9,12,15,3,6,9,12),4,3)
result_test = PPCKO::PPC_KO(data_test,cv,threshold,alpha,n_disc)
result_test$predictions
#you should obtain printed 3.6021  6.6021  9.6021 12.6021
#if yes, everything works


library(fda)
data = t(CanadianWeather$monthlyTemp)

data_used = data[,1:11]
exact_result = data[,12]

start_time_KO <- Sys.time()
pred_KO = PPCKO::PPC_KO(data_used,cv,threshold,alpha,n_disc)
end_time_KO <- Sys.time()

mse_KO = mse(pred_KO$predictions,exact_result)

execution_time_KO <- end_time_KO - start_time_KO
print("Time for KO with C++:")
print(execution_time_KO)
print("MSE for KO with C++:")
print(mse_KO)







start_time_KO_fun <- Sys.time()
rho_KO_fun = PPCdiscretized_KO(data_used)
end_time_KO_fun <- Sys.time()
pred_KO_fun = rho_KO_fun%*%data_used[,11]

mse_KO_fun = mse(pred_KO_fun,exact_result)

execution_time_KO_fun <- end_time_KO_fun - start_time_KO_fun
print("Time for KO with R function")
print(execution_time_KO_fun)
print("MSE for KO with R funciton:")
print(mse_KO_fun)










mse_KO_fun = mse(pred_KO_fun,exact_result)

mse_KO
mse_KO_fun
