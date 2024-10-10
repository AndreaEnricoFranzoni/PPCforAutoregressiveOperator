#include "L2_scalar_product.hpp"

template<typename T>
double
PPC_util::L2_scalar_product<T>::value()
{
  //functional object (a functor, since i have to store the point for which I already have the function evaluation, for interpolating)
  PPC_util::integrand_interp f_integrand{m_grid_evaluations.nodes()};
  
  apsc::NumericalIntegration::Quadrature int_trap(apsc::NumericalIntegration::Trapezoidal{}, m_grid_integration);
  
  return int_trap.apply(f_integrand,m_integrand_evaluations);
}


/*
 * //integrator
 apsc::NumericalIntegration::MonteCarlo integrator;
 apsc::NumericalIntegration::Quadrature mc(integrator, m_grid_integration);
 
 //return the integral value
 return mc.apply(f_integrand,m_integrand_evaluations);
 */