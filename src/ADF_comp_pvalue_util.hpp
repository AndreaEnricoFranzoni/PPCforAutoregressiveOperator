#ifndef ADF_PPC_PVAL_HPP
#define ADF_PPC_PVAL_HPP

#include <vector>
#include <functional>

#include "integrand_interp.hpp"

namespace ADF_util
{

//tables needed to evaluate the p_value
//Critical values for Dickey–Fuller t-distribution.
std::vector<std::vector<double>> table = { {-4.38, -4.15, -4.04, -3.99, -3.98, -3.96},
                                           {-3.95, -3.80, -3.73, -3.69, -3.68, -3.66},
                                           {-3.60, -3.50, -3.45, -3.43, -3.42, -3.41},
                                           {-3.24, -3.18, -3.15, -3.13, -3.13, -3.12},
                                           {-1.14, -1.19, -1.22, -1.23, -1.24, -1.25},
                                           {-0.80, -0.87, -0.90, -0.92, -0.93, -0.94},
                                           {-0.50, -0.58, -0.62, -0.64, -0.65, -0.66},
                                           {-0.15, -0.24, -0.28, -0.31, -0.32, -0.33} };
//time instants of the critical values for Dickey–Fuller t-distribution.
std::vector<double> tableT = {25.0, 50.0, 100.0, 250.0, 500.0, 100000.0};
//honorable cases for pvalue evaluation
std::vector<double> tablep = {0.01, 0.025, 0.05, 0.10, 0.90, 0.95, 0.975, 0.99};


//auxialiary function for the evaluation of the pvalue: interpolate with respect the length of the time series had
std::vector<double>
tableipl(const double &n)
{
  PPC_util::integrand_interp interp_f{tableT};
  
  std::vector<double> tableipl_;
  tableipl_.resize(table.size());
  
  if(n<tableT.front()){std::transform(table.begin(),table.end(),tableipl_.begin(),[](auto el){return el.front();});}
  if(n>tableT.back()) {std::transform(table.begin(),table.end(),tableipl_.begin(),[](auto el){return el.back();});}
  
  std::transform(table.begin(),table.end(),tableipl_.begin(),[&n,&interp_f](auto el){return interp_f(n,el);});
  
  return tableipl_;
}


}

#endif  //ADF_PPC_PVAL_HPP