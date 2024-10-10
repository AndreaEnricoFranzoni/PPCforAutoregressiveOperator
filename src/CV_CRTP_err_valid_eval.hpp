#include "CV_CRTP.hpp"

//error evaluation on a single cv iteration
template<class D, DEF_PARAMS_PPC::cv_strat_type cv_strat, DEF_PARAMS_PPC::cv_err_eval_type err_eval>
double
CV_PPC::CV_KO_PPC_CRTP<D,cv_strat,err_eval>::err_valid_set_eval(const KO_Traits::StoringVector &pred, const KO_Traits::StoringVector &valid, ERR_EVAL_T<DEF_PARAMS_PPC::cv_err_eval_type::MSE>)
const
{
  //using mse between predicted and validation
  return EF_PPC::mse<double>( pred.array() - valid.array() );
}