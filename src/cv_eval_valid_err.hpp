#ifndef CV_EVAL_VALID_ERR_HPP
#define CV_EVAL_VALID_ERR_HPP

#include <algorithm>
#include <numeric>

#include "traits_ko.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif


template<typename T>
double mse(const KO_Traits::StoringVector &diff, int number_threads)
{
  double mse = 0.0;
  int num_comp = diff.size();
  
#ifdef _OPENMP
#pragma omp parallel for shared(diff,num_comp) num_threads(number_threads) reduction(+:mse)
  for(std::size_t i = 0; i < num_comp; ++i)
  {
    mse += std::pow(diff(i),2);
  }
#else
  mse = std::transform_reduce(diff.begin(),
                              diff.end(),
                              0.0,
                              std::plus{},
                              [] (T const &x) {return std::pow(x,2);});
#endif
  
  return mse/static_cast<double>(num_comp);
};

#endif /*CV_EVAL_VALID_ERR_HPP*/