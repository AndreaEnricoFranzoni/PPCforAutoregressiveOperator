// Copyright (c) 2024 Andrea Enrico Franzoni (andreaenrico.franzoni@gmail.com)
//
// This file is part of PPCKO
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of PPCKO and associated documentation files (the PPCKO software), to deal
// PPCKO without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of PPCKO, and to permit persons to whom PPCKO is
// furnished to do so, subject to the following conditions:
//
// PPCKO IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH PPCKO OR THE USE OR OTHER DEALINGS IN
// PPCKO.

#ifndef LIN_REG_PPC_HPP
#define LIN_REG_PPC_HPP

#include <algorithm>
#include <cmath>

#include "traits_ko.hpp"
#include "mlpack/src/mlpack.hpp"
#include "mlpack/src/mlpack/methods/linear_regression/linear_regression.hpp"

/*!
* @file ADF_lr.hpp
* @brief Performing lineat regression. 
* @note The test implementation is as in (https://cran.r-project.org/web/packages/tseries/index.html)
* @author Andrea Enrico Franzoni
*/


/*!
* @class lr_adf
* @brief Class for wrapping Eigen objects into Armadillo ones to perform linear regression
* @tparam LAG_policy indicates if lag orders bigger than one has to be taken into account while computing the test statistic
*/
class lr_adf
{
private:
  /*!Covariates: each column is a statistical unit, each row a covariate: mlpack works as this*/
  KO_Traits::StoringMatrix m_x;         
  /*!Responses*/
  KO_Traits::StoringMatrix m_y;         
  /*!Mse on residuals, between responses and fitted values*/
  double m_mse_res;
  /*!Residual standard error: estimator of sigma_squared (mse/df)=(mse/(statistical units - (covariates(exc. int)+1)))*/                     
  double m_rse;
  /*!Linear regression coefficients*/                         
  KO_Traits::StoringVector m_coeff;     
  /*!Standard errors in the estimate of the coefficients*/
  KO_Traits::StoringVector m_se_coeff;  

  
public:
  /*!
  * @brief Class constructor
  * @param x covariates, stored as column-wise
  * @param y responses
  * @details covariates are transposes in the construction phases, to be coherent with mlpack
  * @note Universal constructor: move semantic used to optimazing handling big size objects
  */
  template<typename STOR_OBJ1, typename STOR_OBJ2>
  lr_adf(STOR_OBJ1&& x, STOR_OBJ2&& y)     
    : m_x{std::forward<STOR_OBJ1>(x.transpose())},m_y{std::forward<STOR_OBJ2>(y)} {}
  

  /*!
  * @brief Standard errors in the estimate of the coefficients: modifying m_se_coeff
  * @param x covariates (stored row-wise)
  */
  inline void std_er_coeff(const arma::mat &x)
  {
    //covariates rank
    auto r = rank(x);

    //full rank: x*x' symmetric and invertible
    if(r == std::min(x.n_rows,x.n_cols))  
    {
      //using Armadillo objects: creating x*x'
      arma::mat inv = arma::inv_sympd(x*trans(x));
      //inverse
      arma::vec diag_inv = arma::diagvec(inv);
      //mapping back to Eigen the standard errors
      m_se_coeff =  Eigen::Map<KO_Traits::StoringVector>(diag_inv.memptr(), diag_inv.n_elem);
      m_se_coeff = m_rse*m_se_coeff;
      //square root
      std::for_each(m_se_coeff.begin(), m_se_coeff.end(), [](auto &el){el=std::sqrt(el);});
    }
    //not full rank: use pseudo inverse
    else                                  
    {
      //using Armadillo objects: creating x*x'
      arma::mat inv = arma::pinv(x*trans(x));
      //pseudo-inverse
      arma::vec diag_inv = arma::diagvec(inv);
      //mapping back to Eigen the standard errors
      m_se_coeff =  Eigen::Map<KO_Traits::StoringVector>(diag_inv.memptr(), diag_inv.n_elem);
      m_se_coeff = m_rse*m_se_coeff;
      //square root
      std::for_each(m_se_coeff.data(), m_se_coeff.data() + m_se_coeff.size(), [](auto &el){el=std::sqrt(el);});
    }
  }
  

  /*!
  * @brief Fitting the linear regression, retaining coeffiecients and standard errors on their estimate
  */
  inline void solve()
  { 
    //convert covariates and responses (Eigen objects) to Armadillo ones in order to use mlpack
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

    //in order to have the standard error of the estimates of the coefficients, it is necessary to include the constant in the design matrix
    //to be done as this in armadillo
    arma::mat covariates_plus_intercept(covariates.n_rows + 1, covariates.n_cols);
    covariates_plus_intercept.row(0).ones();
    covariates_plus_intercept.rows(1, covariates.n_rows) = covariates;
    //standard errors of coefficients estimates
    this->std_er_coeff(covariates_plus_intercept);
  }
  
  /*!
  * @brief Getter for the linear regression coefficients 
  * @return the private m_coeff
  */
  inline KO_Traits::StoringVector coeff() const {return m_coeff;};

  /*!
  * @brief Getter for the linear regression standard errors on the coefficients 
  * @return the private m_se_coeff
  */
  inline KO_Traits::StoringVector se_coeff() const {return m_se_coeff;};
};

#endif //LIN_REG_PPC_HPP