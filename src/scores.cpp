#include "scores.hpp"


void
PPC::scores::evaluating_scores()
{
  m_scores_evaluations.reserve(m_k);
  
  for(const auto &el : m_integrand_evaluations)
  {
    PPC_util::L2_scalar_product<double> scalar_prod((m_grid_evaluations),(m_grid_integration),el.second);
    m_scores_evaluations.emplace_back(scalar_prod.value());
  }
  
}