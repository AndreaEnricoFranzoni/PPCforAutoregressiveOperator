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

#ifndef CV_EVAL_VALID_ERR_HPP
#define CV_EVAL_VALID_ERR_HPP

#include <algorithm>
#include <numeric>

#include "traits_ko.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif


/*!
* @file cv_eval_valid_err.hpp
* @brief Containing the definitions of the different ways of evaluating the validation error
* @author Andrea Enrico Franzoni
*/


/*!
* @brief Template function for the estimate of the L2 norm
* @tparam T type of the data
* @param diff vector containing the differences
* @param number_threads number of threads for OMP
* @return the L2 norm of 'diff'
* @note eventual usage of 'pragma' directive for OMP
*/
template<typename T>
double mse(const KO_Traits::StoringVector &diff, int number_threads)
{
  double mse = 0.0;
  int num_comp = diff.size();
  
  //if OMP: going parallel
#ifdef _OPENMP
#pragma omp parallel for shared(diff,num_comp) num_threads(number_threads) reduction(+:mse)
  for(std::size_t i = 0; i < num_comp; ++i)
  {
    mse += std::pow(diff(i),2);
  }
#else
  // if not OMP: STL algorithms
  mse = std::transform_reduce(diff.begin(),
                              diff.end(),
                              0.0,
                              std::plus{},
                              [] (T const &x) {return std::pow(x,2);});
#endif
  
  return mse/static_cast<double>(num_comp);
};

#endif /*CV_EVAL_VALID_ERR_HPP*/