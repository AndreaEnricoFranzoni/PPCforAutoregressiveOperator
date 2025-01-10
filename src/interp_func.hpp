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

#ifndef INTEGRAND_INTERP_L2_SCALAR_PROD_HPP
#define INTEGRAND_INTERP_L2_SCALAR_PROD_HPP

#include <vector>
#include <functional>

#include "interp1D_util.hpp"

/*!
* @file interp_func.hpp
* @brief Contains a functor to perform unidimensional-univariate function interpolation
* @author Andrea Enrico Franzoni
*/


/*!
* @struct interp_func
* @brief Contains the function to perform unidimensional-univariate function interpolation and the grid for which the function evaluation are already available
*/
struct interp_func
{
  std::vector<double> grid_eval;  ///< Grid for which the function evaluation are already available.

  /*!
  * @brief Perform the interpolation of the function on a given point
  * @param x the abscissa for which the interpolated value is looked for
  * @param evals vector indicating the value of the function in 'grid_eval'
  * @return the interpolated function value
  */
  double operator()(const double &x, const std::vector<double> &evals)
  {
    return apsc::interp1D(grid_eval, evals, x);
  }
};

#endif  //INTEGRAND_INTERP_L2_SCALAR_PROD_HPP