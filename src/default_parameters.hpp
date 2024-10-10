#ifndef KO_DEFAULT_PARAMS_HPP
#define KO_DEFAULT_PARAMS_HPP

#include <Eigen/Dense>

namespace DEF_PARAMS_PPC
{

//CV parameters
//dimension of the grid to make CV on alpha on log scale
constexpr std::size_t dim_grid_alpha = 21;

//lower exponent for the log grid of alphas
constexpr int min_exp_alphas = -10;

//tolerance for CV on k
constexpr double toll_cv_k = 1e-4;


//integration parameters
//left extreme
constexpr double a_interval = 0.0;

//right extreme
constexpr double b_interval = 1.0;

//dimension of the ingtegration gird for using MC
constexpr int dim_grid_int = 250;


//CV train/validation set strategy
enum cv_strat_type
{
  AUGMENTING_WINDOW = 0,  //using augmenting window during CV
};


//CV error evaluation
enum cv_err_eval_type
{
  MSE = 0,  //using mse to evaluate the prediction of the trained model in the validation
};


//removing NaNs
enum MA_type
{
  MR  = 0,      //replacing nans with mean (easily changes the mean of the distribution)
  ZR  = 1,      //replacing nans with 0s (easily changes the sd of the distribution)
};


} //end namespace DEF_PARAMS_PPC

#endif /*KO_DEFAULT_PARAMS_HPP*/