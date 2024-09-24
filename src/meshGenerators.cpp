#include "meshGenerators.hpp"
#include <algorithm>
#include <stdexcept>

namespace Geometry
{
MeshNodes
Uniform::operator()() const
{
  auto const &n = this->M_num_elements;
  auto const &a = this->M_domain.left();
  auto const &b = this->M_domain.right();
  if(n == 0)
    throw std::runtime_error("At least two elements");
  MeshNodes    mesh(n + 1);
  double const h = (b - a) / static_cast<double>(n);
#pragma omp parallel for
  for(auto i = 0u; i < n; ++i)
    mesh[i] = a + h * i;
  mesh[n] = b;
  return mesh;
  
}

} // namespace Geometry
