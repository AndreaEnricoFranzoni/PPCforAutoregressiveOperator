#ifndef CV_PPC_HPP
#define CV_PPC_HPP

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



namespace CV_PPC            //CrossValidation for PPC
{

using singleCV_t = std::function<KO_Traits::StoringVector(KO_Traits::StoringMatrix,double,double,int)>;
using ef_t = std::function<double(KO_Traits::StoringVector)>;


/*!
 * Class for doing CV within the context of PPC KO algo, for only one of the parameters
 */
template<typename T>                         //double if for alpha, int if for k
class CV_KO
{
private:
  KO_Traits::StoringMatrix m_X;
  std::vector<T> m_params;
  std::vector<double> m_errors;           //the error in the position i-th is relative to the parameter in the position i-th in m_params
  std::size_t m_grid_dim;
  T m_param_best;
  singleCV_t m_singleCV;
  ef_t m_ef;
  
  //parameters needed for making KO working
  double m_threshold_ppc;
  double m_alpha;
  int m_k;
  

public:
  CV_KO(KO_Traits::StoringMatrix&& X,
        std::size_t grid_dim,
        double threshold_ppc,
        double alpha,
        int k,
        const singleCV_t & singleCV,
        const ef_t & ef)
    :   
    m_X{std::forward<KO_Traits::StoringMatrix>(X)},
    m_errors(grid_dim,static_cast<double>(0)),
    m_grid_dim(grid_dim),
    m_threshold_ppc(threshold_ppc),
    m_alpha(alpha),
    m_k(k),
    m_singleCV(singleCV),
    m_ef(ef)
  
    {
      m_params.resize(m_grid_dim);
      
      if constexpr(std::is_same<T,int>::value)      //CV su k: try values from 1 to m
      {
        std::iota(m_params.begin(),m_params.end(),static_cast<T>(1));
      }
      
      if constexpr(std::is_same<T,double>::value)   //CV su alpha 
      { 
        std::iota(m_params.begin(),m_params.end(),static_cast<T>(DEF_PARAMS_PPC::min_exp_alphas));
        std::transform(m_params.begin(),m_params.end(),m_params.begin(),[](T el){return(pow(static_cast<T>(10),el));});
      }
    }
  

  /*!
   * Getter for m_X
   */
  inline KO_Traits::StoringMatrix X() const {return m_X;};
  
  /*!
   * Getter for m_param_best
   */
  inline T param_best() const {return m_param_best;};
  
  /*!
   * Getter for m_singleCV
   */
  inline singleCV_t singleCV() const {return m_singleCV;};
  
  /*!
   * Getter for m_ef
   */
  inline ef_t ef() const {return m_ef;};
  
  
  
  /*!
   * Getter for m_errors (only needed for debugging)
   */
  inline std::vector<double> errors() const {return m_errors;};
  
  
  //doing a single CV, for a given param (returns the error evaluated using the error function m_ef for a given parameter and a given training set)
  double single_cv(T param, int dim_train_set) const;
  
  //returns the mean of the m_ef for a given parameter evaluated on the different training set (obtained moving the window, adding every time one time instant, and using the following one as validation set)
  double moving_window_cv(T param, const std::vector<int> & t_i) const;
  
  //find best param according to average of mse of the predictions on the moving window
  void best_param();
};

}   //end namespace CV_PPC


#include "CV_KO_imp.hpp"

#endif  /*CV_PPC_HPP*/