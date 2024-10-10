#include <RcppEigen.h>


#include "default_parameters.hpp"
#include "wrapper_params.hpp"

#include "reading_data.hpp"
#include "PPC_KO.hpp"
#include "KO_factory.hpp"

#include "ADF_lr.hpp"
#include "ADF_test.hpp"
#include "ADF_policies.hpp"


using namespace Rcpp;


//
// [[Rcpp::depends(RcppEigen)]]



//
// [[Rcpp::export]]
Rcpp::List PPC_KO(Rcpp::NumericMatrix X,
                  std::string id_CV = "NoCV",
                  double threshold_ppc = 0.95,
                  double alpha = 0.75,
                  int k = 0,
                  Rcpp::Nullable<std::string> id_rem_nan = R_NilValue
                  )
{ 
  using T = double;       //version for real-values time series
  
  //wrapping parameters
  const DEF_PARAMS_PPC::MA_type id_RN = WRAP_PARAMS_PPC::wrap_id_rem_nans(id_rem_nan);
  
  
  //reading data, handling NANs
  KO_Traits::StoringMatrix x = reading_data<T>(X,id_RN);
  
  
  //checking parameters
  WRAP_PARAMS_PPC::check_threshold_ppc(threshold_ppc);
  WRAP_PARAMS_PPC::check_alpha(alpha);
  WRAP_PARAMS_PPC::check_k(k,x.rows());
  
  
  
  //object solver  !!using std::move() to move the data (better if there are big data)!!
  std::unique_ptr<PPC::PPC_KO_base> ko = KO_Factory::KO_solver(id_CV,std::move(x),threshold_ppc,alpha,k);
  
  //solving using the requested KO algo
  ko->solve();

  //estimate of the predictions
  auto one_step_ahead_pred  = ko->prediction();
  //alpha used
  double alpha_used         = ko->alpha();
  //number of PPC retained
  int n_PPC                 = ko->k();
  //scores along the k PPCs
  auto scores_PPC           = ko->scores();
  //estimate of rho
  auto rho_estimate         = ko->rho();
  
  auto valid_err = ko->ValidErr();
  

  //saving results in a list
  Rcpp::List l;
  l["predictions"]    = one_step_ahead_pred;
  l["alpha"]          = alpha_used;
  l["PPCs_retained"]  = n_PPC;
  l["scores"]         = scores_PPC;
  l["rho_hat"]        = rho_estimate;
  l["valid_errors"]   = valid_err;
  
  
  

  
  return l;
}




//
// [[Rcpp::export]]
Rcpp::List read_data_na(Rcpp::NumericMatrix X)
{ 
  using T = double;       //version for double
  
  DEF_PARAMS_PPC::MA_type id_RN = DEF_PARAMS_PPC::MA_type::MR;
  KO_Traits::StoringMatrix x = reading_data<T>(X, id_RN);
  
  
  //saving results in a list
  Rcpp::List l;
  l["data_read"] = x;

  return l;
}


//
// [[Rcpp::export]]
Rcpp::List test_lr(Rcpp::NumericMatrix X, Rcpp::NumericMatrix Y)
{
  using T = double; 
  
  KO_Traits::StoringMatrix x = reading_data<T>(X, DEF_PARAMS_PPC::MA_type::MR);
  KO_Traits::StoringMatrix y = reading_data<T>(Y, DEF_PARAMS_PPC::MA_type::MR);
  
  ADF_util::lr_adf rg(std::move(x),std::move(y));
  
  rg.solve();
  
  auto coeff = rg.coeff();
  auto se_coeff = rg.se_coeff();
  
  Rcpp::List l;
  l["coeff"] = coeff;
  l["se_coeff"] = se_coeff;
  return l;
}


//
// [[Rcpp::export]]
Rcpp::List test_adf(Rcpp::NumericMatrix X, int lag_order)
{
  using T = double; 
  
  KO_Traits::StoringMatrix x = reading_data<T>(X, DEF_PARAMS_PPC::MA_type::MR);
  
  
  
  PPC_util::adf<ADF_util::CaseLagOrderADF> adf_t(std::move(x),lag_order);
  
  adf_t.test();
  auto pv = adf_t.p_values();
  
  
  Rcpp::List l;
  l["pvalues"] = pv;
  return l;
}



/*
 * m = dim(X.sample)[1]
 p_vals = numeric(m)
 
 for(i in 1:m){
 l = adf.test(X.sample[i,])
 p_vals[i] = l$p.value
 }
 
 c_b = PPCKO::test_adf(X.sample,4)
 
 PPCKO::PPC_KO(X.sample,"CV_k")
 
 
 
 
#test for adf
 X.test = as.matrix(X.sample[1:3,])
 X.test = matrix(data=c(1,2,4,7,11,16,22,29,37,46,56), nrow=2, ncol=length(c(1,2,4,7,11,16,22,29,37,46,56)), byrow = T)
 diff(X.test)
 
 
 x_dii=diff(X.test[3,])
 t=PPCKO::read_data_na(X.test)
 
 
 
 x = X.test[3,]
 y = diff(x = x)
 k_start = 4
 k = k_start + 1
 n = length(y)
 z = embed(y,k)
 z
 yt = z[,1]
 yt
 xt1 = x[k:n]
 xt1
 tt = k:n
 tt
 yt1 <- z[,2:k]
 res <- lm(yt ~ xt1 + 1 + tt + yt1)
 res.sum <- summary(res)
 res.sum
 res.sum$coefficients[2,1]/res.sum$coefficients[2,2] 
 h = PPCKO::test_adf(X.test,k_start)
 l = adf.test(x,alternative = "stationary")
 l$alternative
 
 
 
 
 STAT = res.sum$coefficients[2,1]/res.sum$coefficients[2,2] 
 
 table <- cbind(c(4.38, 4.15, 4.04, 3.99, 3.98, 3.96),
 c(3.95, 3.80, 3.73, 3.69, 3.68, 3.66),
 c(3.60, 3.50, 3.45, 3.43, 3.42, 3.41),
 c(3.24, 3.18, 3.15, 3.13, 3.13, 3.12),
 c(1.14, 1.19, 1.22, 1.23, 1.24, 1.25),
 c(0.80, 0.87, 0.90, 0.92, 0.93, 0.94),
 c(0.50, 0.58, 0.62, 0.64, 0.65, 0.66),
 c(0.15, 0.24, 0.28, 0.31, 0.32, 0.33))
 table <- -table
 tablen <- dim(table)[2]
 
 tableT <- c(25, 50, 100, 250, 500, 100000)
 tablep <- c(0.01, 0.025, 0.05, 0.10, 0.90, 0.95, 0.975, 0.99)
 tableipl <- numeric(tablen)
 for(i in (1:tablen)){tableipl[i] <- approx(tableT, table[, i], n, rule=2)$y}
 
 
 interpol <- approx(tableipl, tablep, STAT, rule=2)$y
 
 
 
 
 
 
#test for regression
 data=cars
 data = data.frame(data[,1],data[,1]^2,data[,1]^3,data[,2])
 colnames(data)=c("speed","speed^2","speed^3","dist")
 n = dim(data)[2]
 
 
 
 lr_R = lm( data[,n] ~ data[,1]+1+data[,2]+data[,3])
 
 res.sum <- summary(lr_R)
 res.sum$coefficients[2,1]
 res.sum
 
 design_matrix = model.matrix(lr_R)
 solve(t(design_matrix)%*%design_matrix)
 
 
 cov = as.matrix(data[,1:(n-1)])
 resp = as.matrix(data[,n])
 lr_C = PPCKO::test_lr(cov,resp)
 
 
 
#Residual Standard error (Like Standard Deviation)
 k=3 #Subtract one to ignore intercept
 SSE=sum(lr_R$residuals**2)
 n=length(lr_R$residuals)
 sqrt(SSE/(n-(1+k)))/sqrt(sum(data[,1]^2))
 sqrt(SSE/(n-(1+k))*solve(t(X) %*% X)[2,2])
 solve(t(design_matrix)%*%design_matrix)
 
 
 
 
 
 
 
 
 
 
 */