#include "CV.hpp"

//error evaluation on a single cv iteration
template< class D, CV_STRAT cv_strat, CV_ERR_EVAL err_eval, K_IMP k_imp, VALID_ERR_RET valid_err_ret >
double
CV_base<D,cv_strat,err_eval,k_imp,valid_err_ret>::err_valid_set_eval(const KO_Traits::StoringVector &pred, const KO_Traits::StoringVector &valid, int number_threads, ERR_EVAL_T<CV_ERR_EVAL::MSE>)
const
{
  //using mse between predicted and validation
  return mse<double>( pred.array() - valid.array(), number_threads );
}