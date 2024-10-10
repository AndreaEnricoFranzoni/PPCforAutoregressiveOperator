#ifndef ADF_PPC_POLICIES_HPP
#define ADF_PPC_POLICIES_HPP

#include <algorithm>
#include <numeric>
#include <array>

#include "KO_Traits.hpp"
#include "ADF_lr.hpp"


namespace ADF_util
{

//returning an array containing the two values needed for computing adf test statistics

//no lag order for computing the statistics
struct CaseNoLagOrderADF
{  //n is the number of time differences
  //x: time serie                         dim (n+1) x 1
  //z: embedded time serie difference     dim n x 1
  std::array<double,2> 
  operator()(const KO_Traits::StoringVector &x, const KO_Traits::StoringMatrix &z) const
  {
    //design matrix
    KO_Traits::StoringMatrix covariates(x.size()-1,2);
    covariates.col(0) = x.head(x.size()-1);
    std::iota(covariates.col(1).begin(),covariates.col(1).end(),static_cast<double>(1.0));
    
    //linear regression
    ADF_util::lr_adf lr(std::move(covariates),std::move(z.col(0)));
    lr.solve();
    
    //returning the values necessary for the test statistics
    std::array<double,2> a = {lr.coeff()(1),lr.se_coeff()(1)};
    return a;
  }
};


//lag order for computing the statistic
struct CaseLagOrderADF
{
  //n is the number of time differences
  //x: time serie                         dim (n+1) x 1
  //z: embedded time serie difference     dim n+1 -(k_used) x k_used
  std::array<double,2>
  operator()(const KO_Traits::StoringVector &x, const KO_Traits::StoringMatrix &z) const
  {
    //dimension
    auto k = z.cols();                    //dimension of the embedding space (k_used)
    auto n_row_units = x.size() - k;      //the actual number of statistical units in this regression
    
    //design matrix
    KO_Traits::StoringMatrix covariates(n_row_units,k+1);
    covariates.col(0) = x.segment(k-1,n_row_units);
    std::iota(covariates.col(1).begin(),covariates.col(1).end(),static_cast<double>(k));
    covariates.rightCols(k-1) = z.rightCols(k-1);
    
    //linear regression
    ADF_util::lr_adf lr(std::move(covariates),std::move(z.col(0)));
    lr.solve();
    
    //returning the values necessary for the test statistics
    std::array<double,2> a = {lr.coeff()(1),lr.se_coeff()(1)};
    return a;
  }
};

}   //end namespace ADF_util

#endif  //ADF_PPC_POLICIES_HPP