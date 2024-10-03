data=cars
data = data.frame(data[,1],data[,1]^2,data[,1]^3,data[,2])
colnames(data)=c("speed","speed^2","speed^3","dist")
n = dim(data)[2]



lr_R = lm( data[,n] ~ data[,1]+1+data[,2]+data[,3])

res.sum <- summary(lr_R)
res.sum$coefficients[2,1]


cov = as.matrix(data[,1:(n-1)])
resp = as.matrix(data[,n])
lr_C = PPCKO::test_lr(cov,resp)
