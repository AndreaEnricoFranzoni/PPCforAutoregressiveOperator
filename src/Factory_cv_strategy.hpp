#ifndef CV_STRAT_FACTORY_HPP
#define CV_STRAT_FACTORY_HPP


#include <memory>
#include <utility>

#include "traits_ko.hpp"
#include "strategy_cv.hpp"



//template <typename Integrator, unsigned int ORDER, unsigned int mydim, unsigned int ndim>
//! A Factory class: A class for the choice of the cross-validation method to use for the regularization parameter alpha and number of PPCs k
template<CV_STRAT cv_strat> 
class Factory_cv_strat
{	
public:
  //! Static method that takes a string as identifier and builds a pointer to the right object for the cross-validation requested
  template<typename... Args>
  static 
  std::unique_ptr<cv_strategy<cv_strat>> 
  cv_strat_obj(Args&&...args)
    {
      //augmenting window
      if constexpr(cv_strat == CV_STRAT::AUGMENTING_WINDOW)
        { return std::make_unique<cv_strategy<cv_strat>>(std::forward<Args>(args)...);}
    }
};

#endif //CV_STRAT_FACTORY_HPP