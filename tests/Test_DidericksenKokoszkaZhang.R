#functions and utilities needed
source("data/far_1_generation/requirements.R")
source("data/far_1_generation/far_1.R")
source("data/far_1_generation/functions.R")
source("data/far_1_generation/prediction_error.R")


id_CV_ko = "CV"           #Ko algorithm

#data: look at "far_1.R" for this parameters (only four kernels, two norm, three errors implemented)
id_kernel <- "gaussian"   #way of generating data 
norm      <- 0.8      #Kernel constant (for the L2 norm of the kernel that has to be <1)
id_noise  <- "1"      #error of the FAR(1) process

proc = feat_far_1_process(id_kernel,norm)
id_kernel   <- proc$kernel
a           <- proc$constant
name_kernel <- proc$name


## Parameters --------------------------------------------------------------------------
n <- 100                  #time instants of the functional temporal series
burnin <- 50              #burnin iterations for FAR(1)
N <- n - burnin           

#grid for the FAR(1)
t.grid <- seq(0,1, length.out=200)
s.grid <- seq(0,1, length.out=200)
grid <- expand.grid(t.grid, s.grid)


## Simulate the data
set.seed(23032000)


## Monte Carlo simulation ------------------------------------------------------------

{ 
  #errors vector
  err.dPPC.KO_en <- numeric(N)     #PF:  Kargin-Onatski PPC
  err.dPPC.KO_rn <- numeric(N) 
  err.dEK_en     <- numeric(N)     #EK:  Kokoszka Estimated Kernel
  err.dEK_rn     <- numeric(N)
  err.dEKI_en    <- numeric(N)     #EKI: Kokoszka Estimated Kernel Improved
  err.dEKI_rn    <- numeric(N)
  err.perf_en    <- numeric(N)     #EX:  exact prediction
  err.perf_rn    <- numeric(N)
  err.mean_en    <- numeric(N)     #MP:  mean prediction
  err.mean_rn    <- numeric(N)
  err.naive_en   <- numeric(N)     #NP:  naive prediction
  err.naive_rn   <- numeric(N)
  
  
  ## 1. Simulate a stationary FAR(1) process according to a specific kernel
  X.sample <- far_1(kernel_id = id_kernel, noise_id = id_noise, n = n, t.grid = t.grid, a = a, burnin = burnin)
  #to predict data with the exact value
  exact_pred = innovation(id_kernel)
  
  ## 2. Center the observations
  X.eval <- t(scale(t(X.sample), center=TRUE, scale=FALSE))
  Xmean.eval <- rowMeans(X.sample)
  
  
  pb <- progress::progress_bar$new(
    format = " MC simulation [:bar] :percent in :elapsed",
    total = N, clear = FALSE, width= 60)
  for(b in 1:N) #b=1
  {
    #train and test set (also are necessary with data non centered)
    X.train <- X.eval[,b:(N-1+b)]
    X.test <- X.eval[,N+b]
    X.train_no_cent = X.sample[,b:(N-1+b)]
    X.test_no_cent = X.sample[,N+b]
    
    
    ## 3. Estimate Psi with different methods
    KO_algo        <- PPCKO::PPC_KO( X = X.train_no_cent, id_CV = id_CV_ko, id_p_for_k = "No", threshold_k = 0.95)
    Psihat.dPPC.KO <- KO_algo$rho_hat   
    Psihat.dEK     <- EKdiscretized(X=X.train,p=3)
    Psihat.dEKI    <- EKimproved(X=X.train,p=3)
    
    
    ## 4. Evaluate the error on the test function
    #PF
    Xhat.dPPC.KO      <- KO_algo$predictions
    err.dPPC.KO_en[b] <- En(X.test_no_cent,Xhat.dPPC.KO,t.grid)
    err.dPPC.KO_rn[b] <- Rn(X.test_no_cent,Xhat.dPPC.KO,t.grid)

    #EK
    Xhat.dEK      <- Psihat.dEK %*% X.train[,N] + Xmean.eval
    err.dEK_en[b] <- En(X.test_no_cent,Xhat.dEK,t.grid)
    err.dEK_rn[b] <- Rn(X.test_no_cent,Xhat.dEK,t.grid)
    
    #EKI
    Xhat.dEKI   <- Psihat.dEKI %*% X.train[,N] + Xmean.eval
    err.dEKI_en[b] <- En(X.test_no_cent,Xhat.dEKI,t.grid)
    err.dEKI_rn[b] <- Rn(X.test_no_cent,Xhat.dEKI,t.grid)

    #EX
    Xhat.perf   <- exact_pred(y=X.sample[,(N-1+b)], t.grid=t.grid, a=a)
    err.perf_en[b] <- En(X.test,Xhat.perf,t.grid)
    err.perf_rn[b] <- Rn(X.test,Xhat.perf,t.grid)
   
    #MP (data centered: the prediction is done with their mean, that is 0)
    err.mean_en[b]  <- En(X.test,rep(0,length(t.grid)),t.grid)
    err.mean_rn[b]  <- Rn(X.test,rep(0,length(t.grid)),t.grid)
    
    #NP (prediction is done using always the last observation)
    err.naive_en[b] <- En(X.test,X.eval[,n],t.grid)
    err.naive_rn[b] <- Rn(X.test,X.eval[,n],t.grid)
    
    pb$tick()
  }
}



###################
## BoxPlot of En ##
###################
err_en <- c(err.naive_en, err.perf_en, err.mean_en, err.dPPC.KO_en, err.dEK_en, err.dEKI_en)
method <- rep(c("naive","perfect", "mean","PPC-KO", "EK", "EKI"), each=N)
En <- data.frame(method, err_en)
method_order<- c("naive", "perfect", "mean", "PPC-KO", "EK", "EKI")
En.box <- En %>% mutate(method=factor(x=method, levels=method_order))

##grouped boxplot
quartz()
pgplot <- ggplot(En.box, aes(x=method, y=err_en, fill=method)) + 
  geom_boxplot() + ggtitle(name_kernel)
pgplot <- pgplot +
  #scale_y_continuous(limits=c(0,0.1)) +
  theme_bw() + 
  labs(x="", y="En", fill = "Prediction method") +
  #labs(x="", y=TeX(r'($\frac{1}{N} \; \sum_{j=1}^N (f_{t+1, j}^b - \hat{f}_{t+1,j}^b)^2$)'), fill="Prediction method") +
  theme(plot.title = element_text(face="bold", hjust=0.5, size=22),
        axis.text.x = element_text(size=22),
        axis.title.x = element_text(size=22),
        axis.text.y = element_text(size=22),
        axis.title.y = element_text(size=22),
        legend.title = element_text(size=22),
        legend.text = element_text(size=22),
        legend.position="bottom",
        legend.direction = "horizontal") +
  guides(fill=guide_legend(nrow=1, byrow=TRUE))
pgplot + 
  theme(legend.position="none")



###################
## BoxPlot of Rn ##
###################
err_rn <- c(err.naive_rn, err.perf_rn, err.mean_rn, err.dPPC.KO_rn, err.dEK_rn, err.dEKI_rn)
method <- rep(c("naive","perfect", "mean","PPC-KO", "EK", "EKI"), each=N)
Rn <- data.frame(method, err_rn)
method_order<- c("naive", "perfect", "mean", "PPC-KO", "EK", "EKI")
Rn.box <- En %>% mutate(method=factor(x=method, levels=method_order))

##grouped boxplot
quartz()
pgplot <- ggplot(Rn.box, aes(x=method, y=err_rn, fill=method)) + 
  geom_boxplot()  + ggtitle(name_kernel)
pgplot <- pgplot +
  #scale_y_continuous(limits=c(0,0.1)) +
  theme_bw() + 
  labs(x="", y="Rn", fill = "Prediction method") +
  #labs(x="", y=TeX(r'($\frac{1}{N} \; \sum_{j=1}^N (f_{t+1, j}^b - \hat{f}_{t+1,j}^b)^2$)'), fill="Prediction method") +
  theme(plot.title = element_text(face="bold", hjust=0.5, size=22),
        axis.text.x = element_text(size=22),
        axis.title.x = element_text(size=22),
        axis.text.y = element_text(size=22),
        axis.title.y = element_text(size=22),
        legend.title = element_text(size=22),
        legend.text = element_text(size=22),
        legend.position="bottom",
        legend.direction = "horizontal") +
  guides(fill=guide_legend(nrow=1, byrow=TRUE))
pgplot + 
  theme(legend.position="none")



#mean and standard error of the Ens
mean(err.dPPC.KO_en)
sqrt(var(err.dPPC.KO_en)/length(err.dPPC.KO_en))
mean(err.dEK_en)
sqrt(var(err.dEK_en)/length(err.dEK_en))
mean(err.dEKI_en)
sqrt(var(err.dEKI_en)/length(err.dEKI_en))
mean(err.perf_en)
sqrt(var(err.perf_en)/length(err.perf_en))
mean(err.mean_en)
sqrt(var(err.mean_en)/length(err.mean_en))
mean(err.naive_en)
sqrt(var(err.naive_en)/length(err.naive_en))



#mean and standard error of the Rns
mean(err.dPPC.KO_rn)
sqrt(var(err.dPPC.KO_rn)/length(err.dPPC.KO_rn))
mean(err.dEK_rn)
sqrt(var(err.dEK_rn)/length(err.dEK_rn))
mean(err.dEKI_rn)
sqrt(var(err.dEKI_rn)/length(err.dEKI_rn))
mean(err.perf_rn)
sqrt(var(err.perf_rn)/length(err.perf_rn))
mean(err.mean_rn)
sqrt(var(err.mean_rn)/length(err.mean_rn))
mean(err.naive_rn)
sqrt(var(err.naive_rn)/length(err.naive_rn))
