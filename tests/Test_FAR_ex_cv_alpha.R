#Kernels
#Norm(rho) = 1/2
#identity kernel
kernel_id = function(t,s)
{
  c = 0.5
  return(c)
}

#sloping plane (t) kernel
kernel_slop_plane_t = function(t,s)
{
  c = 0.5
  return(c*t)
}

#sloping plane (s) kernel
kernel_slop_plane_s = function(t,s)
{
  c = 0.5
  return(c*s)
}

#gaussian kernel
kernel_gauss = function(t,s)
{
  c = 0.683
  res = c*exp(-(t^2)/2)*exp(-(s^2)/2)
  return(res)
}


#generate the data
phi = kernel_gauss    #kernel used
phi = kernel_id
phi = kernel_slop_plane_t
phi = kernel_slop_plane_s


m = 50    #number of temporal series
n = 100   #number of time instants


#raw data
data = matrix(rep(0.0,m*n),nrow=m,ncol=n)

# discretize function space on [0,1]
t = seq(from=0, to=1, length.out=m)

#initial value
data[,1] = sin(2*pi*t)


for (n_t in seq(from=2,to=n,by=1)) 
{
  
  #innovation
  eps_1 = rnorm(1)
  eps_2 = rnorm(1)
  lambda = 0.5
  err = eps_1*sqrt(2)*sin(2*pi*t) + eps_2*sqrt(2*lambda)*cos(2*pi*t)
  
  #autoregressive part
  rec_term = rep(0,m)
  for (i in seq(from=1,to=m,by=1)) 
  {
    rec_term[i] = sum(phi(t[i],t)*data[,n_t-1])*(t[2]-t[1])
  }
  data[,n_t] = rec_term + err
}



quartz()
#change here the directory
setwd("/Users/andreafranzoni/Documents/Politecnico/Magistrale/PPCforAutoregressiveOperator_local")

set.seed(230300)
#leave these parameters like that for now
cv = "CV_alpha"
n_disc = 100
alpha_min = 0.00000000001
alpha_max = 10
alphas = seq(alpha_min,alpha_max,length.out=n_disc)

#threshold = 0.95
threshold = 0.75



#raw data
data_used = data[,1:(n-1)]
exact_result = data[,n]







#c++ based package
#number of PPC taken from reg. covariance
start_time_KO_rc <- Sys.time()
KO_rc = PPCKO::PPC_KO( X = data_used,
                       id_CV = cv,
                       id_p_for_k = "Yes",
                       threshold_k = threshold,
                       n_disc = n_disc,
                       alpha_min = alpha_min,
                       alpha_max = alpha_max
)
pred_KO_rc = KO_rc$predictions
alpha_rc = KO_rc$alpha
k_rc = KO_rc$PPCs_retained
end_time_KO_rc <- Sys.time()
#to here with only one command to measure computational time properly

mse_KO_rc = mse(pred_KO_rc,exact_result)
execution_time_KO_rc <- end_time_KO_rc - start_time_KO_rc


#number of PPC taken from phi
start_time_KO_phi <- Sys.time()
KO_phi = PPCKO::PPC_KO( X = data_used,
                        id_CV = cv,
                        id_p_for_k = "No",
                        threshold_k = threshold,
                        n_disc = n_disc,
                        alpha_min = alpha_min,
                        alpha_max = alpha_max
)
pred_KO_phi = KO_phi$predictions
alpha_phi = KO_phi$alpha
k_phi = KO_phi$PPCs_retained
end_time_KO_phi <- Sys.time()
#to here with only one command to measure computational time properly

mse_KO_phi = mse(pred_KO_phi,exact_result)
execution_time_KO_phi <- end_time_KO_phi - start_time_KO_phi


#R func
start_time_KO_fun <- Sys.time()
alpha_opt_fun = dalpha_select(alphas,data_used)$alpha.opt
rho_KO_fun = PPCdiscretized_KO(data_used,p.thresh=threshold,alpha = alpha_opt_fun)
pred_KO_fun = rho_KO_fun%*%data_used[,n-1]
end_time_KO_fun <- Sys.time()

mse_KO_fun = mse(pred_KO_fun,exact_result)
execution_time_KO_fun <- end_time_KO_fun - start_time_KO_fun




#functional data
time <- 1:n
nbasis = ceiling(n/3)
basis <- create.fourier.basis(rangeval=c(1,n),nbasis=nbasis)
basismat <- eval.basis(time, basis) 
est_coef = lsfit(basismat, t(data), intercept = FALSE)$coef
Xsp0 <- basismat %*% est_coef 
data_f = t(Xsp0)
row.names(data_f) = row.names(data)
colnames(data_f) = colnames(data)

#quartz()
#x11()
#plot.fd(Data2fd(y = t(data),argvals = time,basisobj = basis),titles = row.names(data))

data_used_f = data_f[,1:(n-1)]
exact_result_f = data_f[,n]






#c++ based package
#number of PPC taken from reg. covariance
start_time_KO_rc_f <- Sys.time()
KO_rc_f = PPCKO::PPC_KO( X = data_used_f,
                         id_CV = cv,
                         id_p_for_k = "Yes",
                         threshold_k = threshold,
                         n_disc = n_disc,
                         alpha_min = alpha_min,
                         alpha_max = alpha_max
)
pred_KO_rc_f = KO_rc_f$predictions
alpha_rc_f = KO_rc_f$alpha
k_rc_f = KO_rc_f$PPCs_retained
end_time_KO_rc_f <- Sys.time()
#to here with only one command to measure computational time properly

mse_KO_rc_f = mse(pred_KO_rc_f,exact_result_f)
execution_time_KO_rc_f <- end_time_KO_rc_f - start_time_KO_rc_f



#number of PPC taken from phi
start_time_KO_phi_f <- Sys.time()
KO_phi_f = PPCKO::PPC_KO( X = data_used_f,
                          id_CV = cv,
                          id_p_for_k = "No",
                          threshold_k = threshold,
                          n_disc = n_disc,
                          alpha_min = alpha_min,
                          alpha_max = alpha_max
)
pred_KO_phi_f = KO_phi_f$predictions
alpha_phi_f = KO_phi_f$alpha
k_phi_f = KO_phi_f$PPCs_retained
end_time_KO_phi_f <- Sys.time()
#to here with only one command to measure computational time properly

mse_KO_phi_f = mse(pred_KO_phi_f,exact_result_f)
execution_time_KO_phi_f <- end_time_KO_phi_f - start_time_KO_phi_f



#R func
start_time_KO_fun_f <- Sys.time()
alpha_opt_fun_f = dalpha_select(alphas,data_used_f)$alpha.opt
rho_KO_fun_f = PPCdiscretized_KO(data_used_f,p.thresh=threshold,alpha = alpha_opt_fun_f)
pred_KO_fun_f = rho_KO_fun_f%*%data_used_f[,n-1]
end_time_KO_fun_f <- Sys.time()

mse_KO_fun_f = mse(pred_KO_fun_f,exact_result_f)
execution_time_KO_fun_f <- end_time_KO_fun_f - start_time_KO_fun_f




print("CV alpha comparison")

print("Raw data")
print("Time for KO with C++ taking number of PPC from reg cov:")
print(execution_time_KO_rc)
print("Time for KO with C++ taking number of PPC from phi:")
print(execution_time_KO_phi)
print("Time for KO with R function:")
print(execution_time_KO_fun)
print("MSE for KO with C++ taking number of PPC from reg cov:")
print(mse_KO_rc)
print("MSE for KO with C++ taking number of PPC from phi:")
print(mse_KO_phi)
print("MSE for KO with R function:")
print(mse_KO_fun)
print("alpha opt for KO with C++ taking number of PPC from reg cov:")
print(KO_rc$alpha)
print("alpha opt for KO with C++ taking number of PPC from phi:")
print(KO_phi$alpha)
print("alpha opt for KO with R function:")
print(alpha_opt_fun)
print("k for KO with C++ taking number of PPC from reg cov:")
print(KO_rc$PPCs_retained)
print("k for KO with C++ taking number of PPC from phi:")
print(KO_phi$PPCs_retained)


print("Functional data")
print("Time for KO with C++ taking number of PPC from reg cov:")
print(execution_time_KO_rc_f)
print("Time for KO with C++ taking number of PPC from phi:")
print(execution_time_KO_phi_f)
print("Time for KO with R function:")
print(execution_time_KO_fun_f)
print("MSE for KO with C++ taking number of PPC from reg cov:")
print(mse_KO_rc_f)
print("MSE for KO with C++ taking number of PPC from phi:")
print(mse_KO_phi_f)
print("MSE for KO with R function:")
print(mse_KO_fun_f)
print("alpha opt for KO with C++ taking number of PPC from reg cov:")
print(KO_rc_f$alpha)
print("alpha opt for KO with C++ taking number of PPC from phi:")
print(KO_phi_f$alpha)
print("alpha opt for KO with R function:")
print(alpha_opt_fun_f)
print("k for KO with C++ taking number of PPC from reg cov:")
print(KO_rc_f$PPCs_retained)
print("k for KO with C++ taking number of PPC from phi:")
print(KO_phi_f$PPCs_retained)


quartz()
par(mfrow=c(2,2))
plot(alphas,KO_rc$valid_errors,main = "C++ based KO from reg cov raw",xlab = "Alpha",ylab="Validation errors",type = "l")
plot(alphas,KO_rc$valid_errors,main = "C++ based KO from phi raw",xlab = "Alpha",ylab="Validation errors", type="l")
plot(alphas,KO_rc_f$valid_errors,main = "From reg cov fun",xlab = "Alpha",ylab="Validation errors", type="l")
plot(alphas,KO_rc_f$valid_errors,main = "From phi f",xlab = "Alpha",ylab="Validation errors", type="l")



#errors
mses_ko = numeric(n_disc)
mses_ko_f = numeric(n_disc)
mses_fun = numeric(n_disc)
mses_fun_f = numeric(n_disc)

for (i in 1:n_disc) {
  pred_loop_ko = PPCKO::PPC_KO(data_used,"NoCV","No",threshold,alphas[i],n_disc)
  pred_loop_ko_f = PPCKO::PPC_KO(data_used_f,"NoCV","No",threshold,alphas[i],n_disc)
  pred_loop_fun = PPCdiscretized_KO(data_used,alpha = alphas[i])%*%data_used[,n-1]
  pred_loop_fun_f = PPCdiscretized_KO(data_used_f,alpha = alphas[i])%*%data_used_f[,n-1]
  
  mses_ko[i] = mse(pred_loop_ko$predictions,exact_result)
  mses_ko_f[i] = mse(pred_loop_ko_f$predictions,exact_result_f)
  mses_fun[i] = mse(pred_loop_fun,exact_result)
  mses_fun_f[i] = mse(pred_loop_fun_f,exact_result_f)
}

#raw data
y_low = min(min(min(mses_ko),min(mses_fun)),min(KO_rc$valid_errors))
y_up = max(max(max(mses_ko),max(mses_fun)),max(KO_rc$valid_errors))
#functional data data
y_low_f = min(min(min(mses_ko_f),min(mses_fun_f)),min(KO_rc_f$valid_errors))
y_up_f = max(max(max(mses_ko_f),max(mses_fun_f)),max(KO_rc_f$valid_errors)) 





#run from there to the end all together
quartz()
par(mfrow=c(2,1))
plot(alphas,mses_ko,col='black',cex=0.5,type = "l",main="Raw data",
     xlab="Alpha",ylab="Errors",ylim=c(y_low,y_up))
abline(v=alphas[which(mses_ko==min(mses_ko))],col='black')
points(alphas,KO_rc$valid_errors,cex=0.5,type = "l",col='green')
abline(v=alphas[which(KO_rc$valid_errors==min(KO_rc$valid_errors))],col='green')
points(alphas,mses_fun,cex=0.5,type = "l",col='blue')
abline(v=alphas[which(mses_fun==min(mses_fun))],col='blue')
legend("topright",
       legend = c("Test error KO C++ based",
                  "Validation error KO C++ based",
                  "Test error KO R fun based"),
       fill  = c("black","green","blue")
)


plot(alphas,mses_ko_f,col='black',cex=0.5,type = "l",main="Functional data",
     xlab="Alpha",ylab="Errors",ylim=c(y_low,y_up))
abline(v=alphas[which(mses_ko_f==min(mses_ko_f))],col='black')
points(alphas,KO_rc_f$valid_errors,cex=0.5,type = "l",col='green')
abline(v=alphas[which(KO_rc_f$valid_errors==min(KO_rc_f$valid_errors))],col='green')
points(alphas,mses_fun_f,cex=0.5,type = "l",col='blue')
abline(v=alphas[which(mses_fun_f==min(mses_fun_f))],col='blue')
legend("topright",
       legend = c("Test error KO C++ based",
                  "Validation error KO C++ based",
                  "Test error KO R fun based"),
       fill  = c("black","green","blue")
)