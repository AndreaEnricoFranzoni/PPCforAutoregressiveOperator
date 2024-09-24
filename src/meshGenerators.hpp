#ifndef HH_GENERATOR_HH
#define HH_GENERATOR_HH
#include "domain.hpp"
#include <functional>
#include <stdexcept>
#include <vector>
namespace Geometry
{
using MeshNodes = std::vector<double>;
//! General interface
class OneDMeshGenerator
{
public:
  OneDMeshGenerator(Geometry::Domain1D const &d) : M_domain{d} {}
  virtual MeshNodes operator()() const = 0;
  Domain1D
  getDomain() const
  {
    return M_domain;
  }
  virtual ~OneDMeshGenerator() = default;

protected:
  Geometry::Domain1D M_domain;
};
/*! \defgroup meshers Functors which generates a 1D mesh.
  @{ */
//! Uniform mesh
class Uniform : public OneDMeshGenerator
{
public:
  /*! constructor
@param domain A 1D domain
@param b num_elements Number of elements
  */
  Uniform(Geometry::Domain1D const &domain, unsigned int num_elements)
    : OneDMeshGenerator(domain), M_num_elements(num_elements)
  {}
  //! Call operator
  /*!
    @param meshNodes a mesh of nodes
  */
  MeshNodes operator()() const override;

private:
  std::size_t M_num_elements;
};
/*! @}*/
} // namespace Geometry
#endif
