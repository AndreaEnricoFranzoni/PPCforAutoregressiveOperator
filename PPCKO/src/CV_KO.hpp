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

#include "KO_Traits.hpp"
#include "PPC_KO.hpp"
#include "error_function.hpp"


namespace CV_PPC            //CrossValidation for PPC
{

/*!
 * Class for doing CV within the context of PPC KO algo
 */
class CV_KO
{
private:
  KO_Traits::StoringMatrix m_X;
  std::vector<double> m_alphas;
  std::vector<double> m_errors;           //the error in the position i-th is relative to the parameter in the position i-th in m_alphas
  double m_alpha_best;
  std::size_t m_n_disc;
  static constexpr double alpha_max = 10.0;
  static constexpr double alpha_min = 0.0;
  
public:
  CV_KO(KO_Traits::StoringMatrix&& X, std::size_t n_disc)
    :   
    m_X{std::forward<KO_Traits::StoringMatrix>(X)},
    m_errors(n_disc,0.0),
    m_n_disc(n_disc)
    { 
      //the values of alpha tested
      const double step = (alpha_max-alpha_min)/static_cast<double>(n_disc-1);
      m_alphas.resize(m_n_disc);
      double count = 0;
      std::generate(m_alphas.begin(),m_alphas.end(),[&step,&count](){count++; return (alpha_min+step*(count-1));});
    }
  
  /*!
   * Getter for m_X
   */
  inline KO_Traits::StoringMatrix X() const {return m_X;};
  
  /*!
   * Getter for m_alpha_best (only one needed)
   */
  inline double alpha_best() const {return m_alpha_best;};
  
  
  
  /*!
   * Getter for m_errors (only needed for debugging)
   */
  inline std::vector<double> errors() const {return m_errors;};
  
  
  
  //doing a single CV, for a given alpha (returns the mse for a given alpha and a given training set)
  double single_cv(double alpha, int dim_train_set,const std::function<double(KO_Traits::StoringVector)> &lf) const;
  
  //for doing a single CV for a given alpha, is necessary to move the window, using different train and test set
  //(returns the mean of the mses for a given alpha)
  double moving_window_cv(double alpha, const std::vector<int> & t_i, const std::function<double(KO_Traits::StoringVector)> &lf) const;
  
  //choose best param, evaluating all the alphas
  void best_param();
};

}   //end namespace CV_PPC

#endif  /*CV_PPC_HPP*/