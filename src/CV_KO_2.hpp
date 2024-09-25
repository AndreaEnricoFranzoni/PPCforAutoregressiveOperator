#ifndef CV_PPC2_HPP
#define CV_PPC2_HPP

#include <iostream>
#include <functional>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <concepts>

#include "KO_Traits.hpp"
#include "default_parameters.hpp"



namespace CV_PPC
{
  
using singleCV_on_k_t = std::function<std::pair<int,double>(KO_Traits::StoringMatrix,double,double)>;


class CV_KO_2
{
private:
  KO_Traits::StoringMatrix m_X;
  std::size_t m_grid_dim;
  std::vector<double> m_alphas;            //all the value of alphas
  std::map<double,int> m_best_pairs;       //each alphas with its best k
  std::vector<double> m_errors;            //the error in the position i-th is relative to the pair alpha-k in the map, where alpha is the i-th in m_alphas, and k is the one that gives the best validation error given the alpha
  std::pair<double,int> m_params_best;      //best pair in absolute
  singleCV_on_k_t m_singleCV_on_k;

  //parameters needed for making KO working
  double m_threshold_ppc;
  double m_alpha;
  int m_k;
  
  
public:
  CV_KO_2(KO_Traits::StoringMatrix&& X,
          std::size_t grid_dim,
          double threshold_ppc,
          const singleCV_on_k_t & singleCV_on_k)
      :   
      m_X{std::forward<KO_Traits::StoringMatrix>(X)},
      m_grid_dim(grid_dim),
      m_errors(grid_dim,static_cast<double>(0)),
      m_threshold_ppc(threshold_ppc),
      m_singleCV_on_k(singleCV_on_k)

    {
      m_alphas.resize(m_grid_dim);
      std::iota(m_alphas.begin(),m_alphas.end(),static_cast<double>(DEF_PARAMS_PPC::min_exp_alphas));
      std::transform(m_alphas.begin(),m_alphas.end(),m_alphas.begin(),[](double el){return(pow(static_cast<double>(10),el));});
    }
  
  /*!
   * Getter for m_X
   */
  inline KO_Traits::StoringMatrix X() const {return m_X;};
  
  /*!
   * Getter for m_params_best
   */
  inline std::pair<double,int> params_best() const {return m_params_best;}; //alpha and k
  
  /*!
   * Getter for m_singleCV
   */
  inline singleCV_on_k_t singleCV_on_k() const {return m_singleCV_on_k;};

  
  //given an alpha, does CV on k: returns the best k given that alpha and its validation error
  double single_cv_brute_force(double alpha);

  //find best pair alpha,k given their validation error on the moving window
  void best_params();
};
  
} //end namespace CV_PPC

#endif  /*CV_PPC2_HPP*/