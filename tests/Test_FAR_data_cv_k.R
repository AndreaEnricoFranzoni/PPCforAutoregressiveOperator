library(fda)
data <- t(CanadianWeather$dailyAv[,,1])     #temperature daily
data = t(CanadianWeather$dailyAv[,,2])      #precipitation daily
data = t(CanadianWeather$monthlyTemp)       #temperature monthly
data = t(CanadianWeather$monthlyPrecip)     #precipitation monthly

m = dim(data)[1]      #number of temporal series
n = dim(data)[2]      #number of time instants

quartz()
set.seed(230300)
#leave these parameters like that for now
cv = "CV_k"
alpha = 0.75
k_s = 1:m

threshold=0.95



#raw data
data_used = data[,1:(n-1)]
exact_result = data[,n]






#c++ based package
#number of PPC taken from reg. covariance
start_time_KO_rc <- Sys.time()
KO_rc = PPCKO::PPC_KO( X = data_used,
                       id_CV = cv,
                       id_p_for_k = "No",
                       threshold_k = threshold,
                       id_p_imposed = "Yes"
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
                        id_p_imposed = "No"
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
k_opt_fun = dK_select(data_used,alpha)$nppc.opt
rho_KO_fun = PPCdiscretized_KO(data_used,p.thresh=threshold,alpha = alpha,K=k_opt_fun)
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
                         id_p_for_k = "No",
                         threshold_k = threshold,
                         id_p_imposed = "Yes"
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
                          id_p_imposed = "No"
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
k_opt_fun_f = dK_select(data_used_f,alpha)$nppc.opt
rho_KO_fun_f = PPCdiscretized_KO(data_used_f,p.thresh=threshold,alpha = alpha,K=k_opt_fun_f)
pred_KO_fun_f = rho_KO_fun_f%*%data_used_f[,n-1]
end_time_KO_fun_f <- Sys.time()

mse_KO_fun_f = mse(pred_KO_fun_f,exact_result_f)
execution_time_KO_fun_f <- end_time_KO_fun_f - start_time_KO_fun_f



print("CV k comparison")

print("Raw data")
print("Time for KO with C++ imposing k from reg cov:")
print(execution_time_KO_rc)
print("Time for KO with C++ imposing k from phi:")
print(execution_time_KO_phi)
print("Time for KO with R function:")
print(execution_time_KO_fun)
print("MSE for KO with C++ imposing k from reg cov:")
print(mse_KO_rc)
print("MSE for KO with C++ imposing k from phi:")
print(mse_KO_phi)
print("MSE for KO with R function:")
print(mse_KO_fun)
print("k for KO with C++ imposing k from reg cov:")
print(KO_rc$PPCs_retained)
print("k for KO with C++ imposing k from phi:")
print(KO_phi$PPCs_retained)
print("k for KO withR func:")
print(k_opt_fun)


print("Functional data")
print("Time for KO with C++ imposing k from reg cov:")
print(execution_time_KO_rc_f)
print("Time for KO with C++ imposing k from phi:")
print(execution_time_KO_phi_f)
print("Time for KO with R function:")
print(execution_time_KO_fun_f)
print("MSE for KO with C++ imposing k from reg cov:")
print(mse_KO_rc_f)
print("MSE for KO with C++ imposing k from phi:")
print(mse_KO_phi_f)
print("MSE for KO with R function:")
print(mse_KO_fun_f)
print("k for KO with C++ imposing k from reg cov:")
print(KO_rc_f$PPCs_retained)
print("k for KO with C++ imposing k from phi:")
print(KO_phi_f$PPCs_retained)
print("k for KO withR func:")
print(k_opt_fun_f)


quartz()
par(mfrow=c(2,2))
plot(k_s,KO_rc$valid_errors,main = "C++ based KO from reg cov raw",xlab = "k",ylab="Validation errors",type = "l")
plot(k_s,KO_phi$valid_errors,main = "C++ based KO from phi raw",xlab = "k",ylab="Validation errors", type="l")
plot(k_s,KO_rc_f$valid_errors,main = "From reg cov fun",xlab = "k",ylab="Validation errors", type="l")
plot(k_s,KO_phi_f$valid_errors,main = "From phi f",xlab = "k",ylab="Validation errors", type="l")





#errors
mses_ko = numeric(length = length(k_s))
mses_ko_f = numeric(length = length(k_s))
mses_fun = numeric(length = length(k_s))
mses_fun_f = numeric(length = length(k_s))

for (i in 1:length(k_s)) {
  pred_loop_ko = PPCKO::PPC_KO(data_used,"NoCV","No",threshold,alpha = alpha,k=k_s[i])
  pred_loop_ko_f = PPCKO::PPC_KO(data_used_f,"NoCV","No",threshold,alpha = alpha,k=k_s[i])
  pred_loop_fun = PPCdiscretized_KO(data_used,alpha = alpha,K=k_s[i])%*%data_used[,n-1]
  pred_loop_fun_f = PPCdiscretized_KO(data_used_f,alpha = alpha,K=k_s[i])%*%data_used_f[,n-1]
  
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
plot(k_s,mses_ko,col='black',cex=0.5,type = "l",main="Raw data",
     xlab="k",ylab="Errors",ylim=c(y_low,y_up))
abline(v=k_s[which(mses_ko==min(mses_ko))],col='black')
points(k_s,KO_rc$valid_errors,cex=0.5,type = "l",col='green')
abline(v=k_s[which(KO_rc$valid_errors==min(KO_rc$valid_errors))],col='green')
points(k_s,mses_fun,cex=0.5,type = "l",col='blue')
abline(v=k_s[which(mses_fun==min(mses_fun))],col='blue')
legend("topright",
       legend = c("Test error KO C++ based",
                  "Validation error KO C++ based",
                  "Test error KO R fun based"),
       fill  = c("black","green","blue")
)


plot(k_s,mses_ko_f,col='black',cex=0.5,type = "l",main="Functional data",
     xlab="Alpha",ylab="Errors",ylim=c(y_low,y_up))
abline(v=k_s[which(mses_ko_f==min(mses_ko_f))],col='black')
points(k_s,KO_rc_f$valid_errors,cex=0.5,type = "l",col='green')
abline(v=k_s[which(KO_rc_f$valid_errors==min(KO_rc_f$valid_errors))],col='green')
points(k_s,mses_fun_f,cex=0.5,type = "l",col='blue')
abline(v=k_s[which(mses_fun_f==min(mses_fun_f))],col='blue')
legend("topright",
       legend = c("Test error KO C++ based",
                  "Validation error KO C++ based",
                  "Test error KO R fun based"),
       fill  = c("black","green","blue")
)