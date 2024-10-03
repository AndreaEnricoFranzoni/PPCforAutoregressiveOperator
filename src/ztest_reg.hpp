#ifndef TEST_REG_HPP
#define TEST_REG_HPP


#include "KO_Traits.hpp"
#include "mlpack/src/mlpack.hpp"
#include "mlpack/src/mlpack/methods/linear_regression/linear_regression.hpp"

#include <iostream>


class reg_test
{
private:
  KO_Traits::StoringMatrix m_x;         //covariate salvate per colonne: ogni colonne è un'unità statistica, ogni riga una covariata perchè mlpack funziona con design matrix invertita
  KO_Traits::StoringMatrix m_y;
  
  KO_Traits::StoringVector m_coeff;
  KO_Traits::StoringVector m_res;
  
  
public:
  reg_test(KO_Traits::StoringMatrix&& x,
           KO_Traits::StoringMatrix&& y
  )     //pass the covariates transpose since it is how mlpack works
    : m_x{std::forward<KO_Traits::StoringMatrix>(x.transpose())},m_y{std::forward<KO_Traits::StoringMatrix>(y)} {}
  
  inline void solve()
  {
    arma::mat covariates(m_x.data(), m_x.rows(), m_x.cols(), false, true);
    arma::rowvec responses(m_y.data(), m_y.rows(), false, true);
    
    mlpack::LinearRegression lr(covariates,responses);
    
    m_coeff = Eigen::Map<KO_Traits::StoringVector>(lr.Parameters().memptr(),
                                                   lr.Parameters().n_rows,
                                                   lr.Parameters().n_cols);
  }
  
  inline KO_Traits::StoringVector coeff() const {return m_coeff;};
  
};


#endif //TEST_REG_HPP