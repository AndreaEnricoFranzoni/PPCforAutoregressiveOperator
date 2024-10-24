#ifndef CV_EVAL_VALID_ERR_HPP
#define CV_EVAL_VALID_ERR_HPP

#include <algorithm>
#include <numeric>

#include "traits_ko.hpp"



template<typename T>
double mse(const KO_Traits::StoringVector &diff)
{
    double mse = std::transform_reduce(
    diff.begin(),
    diff.end(),
    0.0,
    std::plus{},
    [] (T const &x) {return std::pow(x,2);}
  );
  
  return mse/static_cast<double>(diff.size());
};

#endif /*CV_EVAL_VALID_ERR_HPP*/