#include "CV_KO_2.hpp"

double
CV_PPC::CV_KO_2::single_cv_brute_force(double alpha)
{   
  //given alpha, performs CV on k: gives back the best k and its validation error
  std::pair<int,double> p = this->singleCV_on_k()(this->X(),m_threshold_ppc,alpha);
  
  //storing, for a given alpha, the best pair alpha-k
  m_best_pairs.insert(std::make_pair(alpha,p.first));
  
  //return the validation error for, for a given alpha, the best pair alpha-k
  return p.second;
}


void
CV_PPC::CV_KO_2::best_params()
{ 
 
  //for all alphas: do a CV on k, retain the best k for the given alpha and its validation error (mean of mse of predictions on moving windows)
  std::transform(m_alphas.cbegin(),
                 m_alphas.cend(),
                 m_errors.begin(),
                 [this](double const &alpha_i){return this->single_cv_brute_force(alpha_i);});
  
  //best pair: the one associated to the smaller paramter
  double best_alpha = m_alphas[std::distance(m_errors.begin(),std::min_element(m_errors.begin(),m_errors.end()))];
  m_params_best = std::make_pair(best_alpha,m_best_pairs.find(best_alpha)->second);
  
  //free memory
  m_alphas.clear();
  m_best_pairs.clear();
  m_errors.clear();
}