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

set.seed(230300)
#leave these parameters like that for now
cv = "CV"
nocv = "NoCV"
threshold = 0.75
alpha = 0.75
n_disc = 1000
alpha_max = 10

quartz()
#x11()
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



#raw data
data_used = data[,1:(n-1)]
exact_result = data[,n]


#c++ based package
#number of PPC taken from reg. covariance
start_time_KO_rc <- Sys.time()
KO_rc = PPCKO::PPC_KO(data_used,nocv,"Yes",threshold,alpha,n_disc)
pred_KO_rc = KO_rc$predictions
k_rc = KO_rc$PPCs_retained
end_time_KO_rc <- Sys.time()
#to here with only one command to measure computational time properly

mse_KO_rc = mse(pred_KO_rc,exact_result)
execution_time_KO_rc <- end_time_KO_rc - start_time_KO_rc


#number of PPC taken from phi
start_time_KO_phi <- Sys.time()
KO_phi = PPCKO::PPC_KO(data_used,nocv,"No",threshold,alpha,n_disc)
pred_KO_phi = KO_phi$predictions
k_phi = KO_phi$PPCs_retained
end_time_KO_phi <- Sys.time()
#to here with only one command to measure computational time properly

mse_KO_phi = mse(pred_KO_phi,exact_result)
execution_time_KO_phi <- end_time_KO_phi - start_time_KO_phi


#R func
start_time_KO_fun <- Sys.time()
rho_KO_fun = PPCdiscretized_KO(data_used,p.thresh=threshold,alpha = alpha)
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

quartz()
#x11()
plot.fd(Data2fd(y = t(data),argvals = time,basisobj = basis),titles = row.names(data))

data_used_f = data_f[,1:(n-1)]
exact_result_f = data_f[,n]


#c++ based package
#number of PPC taken from reg. covariance
start_time_KO_rc_f <- Sys.time()
KO_rc_f = PPCKO::PPC_KO(data_used_f,nocv,"Yes",threshold,alpha,n_disc)
pred_KO_rc_f = KO_rc_f$predictions
k_rc_f = KO_rc_f$PPCs_retained
end_time_KO_rc_f <- Sys.time()
#to here with only one command to measure computational time properly

mse_KO_rc_f = mse(pred_KO_rc_f,exact_result_f)
execution_time_KO_rc_f <- end_time_KO_rc_f - start_time_KO_rc_f


#number of PPC taken from phi
start_time_KO_phi_f <- Sys.time()
KO_phi_f = PPCKO::PPC_KO(data_used_f,nocv,"No",threshold,alpha,n_disc)
pred_KO_phi_f = KO_phi_f$predictions
k_phi_f = KO_phi_f$PPCs_retained
end_time_KO_phi_f <- Sys.time()
#to here with only one command to measure computational time properly

mse_KO_phi_f = mse(pred_KO_phi_f,exact_result_f)
execution_time_KO_phi_f <- end_time_KO_phi_f - start_time_KO_phi_f


#R func
start_time_KO_fun_f <- Sys.time()
rho_KO_fun_f = PPCdiscretized_KO(data_used_f,p.thresh = threshold,alpha = alpha)
pred_KO_fun_f = rho_KO_fun_f%*%data_used[,n-1]
end_time_KO_fun_f <- Sys.time()

mse_KO_fun_f = mse(pred_KO_fun_f,exact_result_f)
execution_time_KO_fun_f <- end_time_KO_fun_f - start_time_KO_fun_f


print("NO CV comparison")

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
print("k for KO with C++ taking number of PPC from reg cov:")
print(KO_rc_f$PPCs_retained)
print("k for KO with C++ taking number of PPC from phi:")
print(KO_phi_f$PPCs_retained)
