set.seed(230300)
#leave these parameters like that for now
cv = "CV"
nocv = "NoCV"
threshold = 0.95
alpha = 0.75
n_disc = 1000
alpha_max = 10

#change here the directory
setwd("/Users/andreafranzoni/Documents/Politecnico/Magistrale/PACS/PPCforAutoregressiveOperator")

#building the dataset
n_ts = 20
n_ti = 30
last_tim = 10
times = seq(0,last_tim,length.out=n_ti)
data = matrix(nrow=n_ts,ncol=n_ti)
beta = rnorm(n_ts,10,12)
v = rnorm(n_ts,40,20)

for (i in 1:20) {
  for (j in 1:30) {
    data[i,j] = beta[i]*sin(times[j]) + v[i] + rnorm(1,0,0.01)
  }
}



m = dim(data)[1]
n = dim(data)[2]

data_used = data[,1:(n-1)]
exact_result = data[,n]




###################
###   CV TEST   ###
###################
#test with c++ based package
#run from here
start_time_KO_cv <- Sys.time()
KO_cv = PPCKO::PPC_KO(data_used,cv,threshold,alpha,n_disc)
pred_KO_cv = KO_cv$predictions
alpha_opt_ko_cv = KO_cv$alpha
end_time_KO_cv <- Sys.time()
#to here with only one command to measure computational time properly

mse_KO_cv = mse(pred_KO_cv,exact_result)
execution_time_KO_cv <- end_time_KO_cv - start_time_KO_cv


#test with R function
#run from here
start_time_KO_fun_cv <- Sys.time()
alpha.vec = seq(1e-7, alpha_max, length.out = n_disc)    #to avoid singular matrices inversion
alpha_opt_fun_cv = dalpha_select(alpha.vec,data_used)$alpha.opt
rho_KO_fun_cv = PPCdiscretized_KO(data_used,alpha = alpha_opt_fun_cv)
pred_KO_fun_cv = rho_KO_fun_cv%*%data_used[,11]
end_time_KO_fun_cv <- Sys.time()
#to here with only one command to measure computational time properly

mse_KO_fun_cv = mse(pred_KO_fun_cv,exact_result)
execution_time_KO_fun_cv <- end_time_KO_fun_cv - start_time_KO_fun_cv


print("CV comparison")
print("Time for KO with cv with C++:")
print(execution_time_KO_cv)
print("Time for KO with cv with R function:")
print(execution_time_KO_fun_cv)
print("MSE for KO with cv with C++:")
print(mse_KO_cv)
print("MSE for KO with cv with R function:")
print(mse_KO_fun_cv)




######################
###   NO CV TEST   ###
######################
#test with c++ based package
#run from here
start_time_KO <- Sys.time()
KO = PPCKO::PPC_KO(data_used,nocv,threshold,alpha,n_disc)
pred_KO = KO$predictions
alpha_opt_ko = KO$alpha
end_time_KO <- Sys.time()
#to here with only one command to measure computational time properly

mse_KO = mse(pred_KO,exact_result)
execution_time_KO <- end_time_KO - start_time_KO


#test with R function
#run from here
start_time_KO_fun <- Sys.time()
rho_KO_fun = PPCdiscretized_KO(data_used,alpha = alpha)
pred_KO_fun = rho_KO_fun%*%data_used[,11]
end_time_KO_fun <- Sys.time()
#to here with only one command to measure computational time properly

mse_KO_fun = mse(pred_KO_fun,exact_result)
execution_time_KO_fun <- end_time_KO_fun - start_time_KO_fun


print("NO CV comparison")
print("Time for KO without cv with C++:")
print(execution_time_KO)
print("Time for KO without cv with R function:")
print(execution_time_KO_fun)
print("MSE for KO without cv with C++:")
print(mse_KO)
print("MSE for KO without cv with R function:")
print(mse_KO_fun)





#nb
alpha.vec = seq(0.001,alpha_max,length.out = n_disc )
mses_ko = numeric(n_disc)
mses_fun = numeric(n_disc)

for (i in 1:n_disc) {
  pred_loop_ko = PPCKO::PPC_KO(data_used,"NoCV",threshold,alpha = alpha.vec[i],n_disc=n_disc)
  pred_loop_fun = PPCdiscretized_KO(data_used,alpha = alpha.vec[i])%*%data_used[,11]
  mses_ko[i] = mse(pred_loop_ko$predictions,exact_result)
  mses_fun[i] = mse(pred_loop_fun,exact_result)
}


y_low = min(min(min(mses_ko),min(mses_fun)),min(KO_cv$valid_errors))
y_up = max(max(max(mses_ko),max(mses_fun)),max(KO_cv$valid_errors)) 
#run from there to the end all together
quartz()
plot(alpha.vec,mses_ko,col='black',cex=0.5,
     xlab="Alpha",ylab="Errors",ylim=c(y_low,y_up))
abline(v=alpha.vec[which(mses_ko==min(mses_ko))],col='black')
points(alpha.vec,KO_cv$valid_errors,cex=0.5,col='green')
abline(v=alpha.vec[which(KO_cv$valid_errors==min(KO_cv$valid_errors))],col='green')
points(alpha.vec,mses_fun,cex=0.5,col='blue')
abline(v=alpha.vec[which(mses_fun==min(mses_fun))],col='blue')
legend("topright",
       legend = c("Test error KO C++ based",
                  "Validation error KO C++ based",
                  "Test error KO R fun based"),
       fill  = c("black","green","blue")
)



for (i in 1:m) {
  quartz()
  plot(times[1:n-1],data_used[i,],xlim=c(0,last_tim), main = i)
  points(last_tim,exact_result[i],col='red')
  points(last_tim,pred_KO_cv[i],col='blue')
}

