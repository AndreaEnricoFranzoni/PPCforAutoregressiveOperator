#ifndef INTEGRAND_INTERP_L2_SCALAR_PROD_HPP
#define INTEGRAND_INTERP_L2_SCALAR_PROD_HPP

#include <vector>
#include <functional>

#include "interp1D_util.hpp"


//the integrand is an interpolation function: since I need the grid for the interpolation, I relay on a functor
struct interp_func
{
  std::vector<double> grid_eval;
  double operator()(const double &x, const std::vector<double> &evals)
  {
    return apsc::interp1D(grid_eval, evals, x);
  }
};

#endif  //INTEGRAND_INTERP_L2_SCALAR_PROD_HPP