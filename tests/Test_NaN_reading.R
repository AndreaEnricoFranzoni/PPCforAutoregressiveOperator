set.seed(230300)
#leave these parameters like that for now
cv = "CV"
nocv = "NoCV"
threshold = 0.95
alpha = 0.75
n_disc = 1000
alpha_max = 10

#change here the directory
setwd("/Users/andreafranzoni/Documents/Politecnico/Magistrale/PPCforAutoregressiveOperator_local")

#building the dataset
n_ts = 20
n_ti = 30
last_tim = 10
times = seq(0,last_tim,length.out=n_ti)
data = matrix(nrow=n_ts,ncol=n_ti)
beta = rnorm(n_ts,10,12)
v = rnorm(n_ts,40,20)

missing_rows = rbinom(n_ts,1,0.1)
missing_cols = rbinom(n_ti,1,0.1)


for (i in 1:20) {
  for (j in 1:30) {
    
    #data[i,j] = beta[i]*sin(times[j]) + v[i] + rnorm(1,0,0.01)
    
    if(missing_rows[i]==0 && missing_cols[j]==0)
    {
      data[i,j] = beta[i]*sin(times[j]) + v[i] + rnorm(1,0,0.01)
    }
    else
    {
      data[i,j] = NaN
    }
  }
}


m = dim(data)[1]
n = dim(data)[2]


data_used = data[,1:(n-1)]
exact_result = data[,n]
n_nan = length(which(is.na(data_used)))
n_nan


res = PPCKO::read_data_na(data_used)
L = res$data_read
n_nan_L = length(which(is.na(L)))
n_nan_L
test_res = PPCKO::PPC_KO(data_used,cv,threshold,alpha,n_disc,id_rem_nan = "ZR")
