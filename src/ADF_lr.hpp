#ifndef LIN_REG_PPC_HPP
#define LIN_REG_PPC_HPP

#include <algorithm>
#include <cmath>

#include "KO_Traits.hpp"
#include "mlpack/src/mlpack.hpp"
#include "mlpack/src/mlpack/methods/linear_regression/linear_regression.hpp"


namespace ADF_util
{

//class for doing linear regression for adf and taking the right coefficients needed for computing its statistic
class lr_adf
{
private:
  KO_Traits::StoringMatrix m_x;         //covariate salvate per colonne: ogni colonne è un'unità statistica, ogni riga una covariata perchè mlpack funziona con design matrix invertita
  KO_Traits::StoringMatrix m_y;         //responses
  
  double m_mse_res;                     //mse on residuals: between responses and fitted values
  double m_rse;                         //residual standard error: estimator of sigma_squared (mse/df)=(mse/(statistical units - (covariates(exc. int)+1)))
  KO_Traits::StoringVector m_coeff;     //coefficient of linear regression
  KO_Traits::StoringVector m_se_coeff;  //standard error in the estimate of the coefficients

  
public:
  lr_adf(KO_Traits::StoringMatrix&& x, KO_Traits::StoringMatrix&& y)     
    : m_x{std::forward<KO_Traits::StoringMatrix>(x.transpose())},m_y{std::forward<KO_Traits::StoringMatrix>(y)} {}
  //pass the covariates and transpose them since it is how mlpack works
  
  
  //standard errors in the estimate of the coefficients
  inline void std_er_coeff(const arma::mat &x)
  {
    auto r = rank(x);
    if(r == std::min(x.n_rows,x.n_cols))  //full rank: x*x' symmetric and invertible
    {
      arma::mat inv = arma::inv_sympd(x*trans(x));
      arma::vec diag_inv = arma::diagvec(inv);
      m_se_coeff =  Eigen::Map<KO_Traits::StoringVector>(diag_inv.memptr(), diag_inv.n_elem);
      m_se_coeff = m_rse*m_se_coeff;
      std::for_each(m_se_coeff.begin(), m_se_coeff.end(), [](auto &el){el=std::sqrt(el);});
    }
    else                                  //not full rank: use pseudo inverse
    {
      arma::mat inv = arma::pinv(x*trans(x));
      arma::vec diag_inv = arma::diagvec(inv);
      m_se_coeff =  Eigen::Map<KO_Traits::StoringVector>(diag_inv.memptr(), diag_inv.n_elem);
      m_se_coeff = m_rse*m_se_coeff;
      std::for_each(m_se_coeff.data(), m_se_coeff.data() + m_se_coeff.size(), [](auto &el){el=std::sqrt(el);});
    }
  }
  
  
  //solving: coefficients and ste on their estimate
  inline void solve()
  { 
    //convert the inputs to arma objects in order to use mlpack
    arma::mat covariates(m_x.data(), m_x.rows(), m_x.cols(), false, true);
    arma::rowvec responses(m_y.data(), m_y.rows(), false, true);
    
    //train the model of linear regression (including intercept)
    mlpack::LinearRegression lr(covariates,responses);
    
    //for adf: needed estaimate of the coefficients and standard errors of these estimates
    //estimate of the coefficients
    m_coeff = Eigen::Map<KO_Traits::StoringVector>(lr.Parameters().memptr(),
                                                   lr.Parameters().n_rows,
                                                   lr.Parameters().n_cols);
    
    //mse between responses and fitted values
    m_mse_res = lr.ComputeError(covariates, responses);
    
    //rse: mse on residulas divided by df (statistical units - covariates (excluding intercept))
    //estimate of sigma_2
    m_rse = m_mse_res*(static_cast<double>(covariates.n_cols)/static_cast<double>(covariates.n_cols - (covariates.n_rows+static_cast<std::size_t>(1)))); 

    //in order to have the standard error of the estimates of the coefficients, I have to include the constant in the design matrix
    //to be done as this in armadillo
    arma::mat covariates_plus_intercept(covariates.n_rows + 1, covariates.n_cols);
    covariates_plus_intercept.row(0).ones();
    covariates_plus_intercept.rows(1, covariates.n_rows) = covariates;
    //standard errors of coefficients estimates
    this->std_er_coeff(covariates_plus_intercept);
    
  }
  
  inline KO_Traits::StoringVector coeff() const {return m_coeff;};
  inline KO_Traits::StoringVector se_coeff() const {return m_se_coeff;};
};

}   //end namespace ADF_util

#endif //LIN_REG_PPC_HPP