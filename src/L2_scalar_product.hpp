#ifndef L2_SCALAR_PRODUCT_HPP
#define L2_SCALAR_PRODUCT_HPP

#include <vector>
#include <utility>

#include "mesh.hpp"
#include "integrand_interp.hpp"
#include "Adams_rule.hpp"
#include "numerical_integration.hpp"

namespace PPC_util
{
  
template<typename T>
class L2_scalar_product
{
  
  
private:
Geometry::Mesh1D m_grid_evaluations;            //values in the domain for which I already have the function value
Geometry::Mesh1D m_grid_integration;            //grid for doing the scalar product (domain information are into the grid)
std::vector<double> m_integrand_evaluations;    //evaluations of the integrand (discrete product between f and g)

  
public:
  template<typename MESH>
  L2_scalar_product(MESH &&grid_evaluations,
                    MESH &&grid_integration,
                    const std::vector<T> &integrand_evaluations
                    )
    : 
    m_grid_evaluations(std::forward<MESH>(grid_evaluations)),
    m_grid_integration(std::forward<MESH>(grid_integration)),               
    m_integrand_evaluations{integrand_evaluations}
  {}
  
  //return the value of the scalar product
  double value();
  
};

}   //end namespace PPC_util

#include "L2_scalar_product_imp.hpp"

#endif  //L2_SCALAR_PRODUCT_HPP