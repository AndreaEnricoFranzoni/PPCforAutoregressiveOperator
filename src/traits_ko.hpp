#ifndef KO_TRAITS_HPP
#define KO_TRAITS_HPP

#include <Eigen/Dense>
#include <tuple>
#include <vector>
#include <variant>
#include <type_traits>
#include <cmath>


// Types for doing algebric operations
struct KO_Traits
{
public:
  
  // Data are stored in dynamic matrices (sizes are known only at compile time, and easily very larges)
  // of doubles: the only conversion from R types is Numeric -> double. some precision is lost. But double are
  // used since time series are of real numbers.
  using StoringMatrix = Eigen::MatrixXd;
  
  using StoringVector = Eigen::VectorXd;    //col-wise
  
  using StoringArray  = Eigen::ArrayXd;
  
};


//domain dimension
enum DOM_DIM
{
  uni_dim  = 0,      //1D domain
  bi_dim   = 1,      //2D domain
};


//if k is imposed 
enum K_IMP
{
  NO  = 0,     //k is not passed as parameter, has to be found using cumulative explanatory power
  YES = 1,     //k is already known 
};


//if validation error has to be stored and returned (memory saving)
enum VALID_ERR_RET
{
  NO_err   = 0,     //validation errors are not stored and not returned
  YES_err  = 1,     //validation errors are stored and returned
};


//CV train/validation set strategy
enum CV_STRAT
{
  AUGMENTING_WINDOW = 0,  //using augmenting window during CV
};


//CV error evaluation
enum CV_ERR_EVAL
{
  MSE = 0,  //using mse to evaluate the prediction of the trained model in the validation
};




// Results

// Depending on the CV (if on one or two params)
using valid_err_cv_1_t = std::vector<double>;
using valid_err_cv_2_t = std::vector<std::vector<double>>;
using valid_err_variant = std::variant<valid_err_cv_1_t,valid_err_cv_2_t>;

// if errors are saved or not
using results_err_t = std::tuple<KO_Traits::StoringVector, double, int, std::vector<double>, KO_Traits::StoringMatrix, valid_err_variant>; 
using results_no_err_t = std::tuple<KO_Traits::StoringVector, double, int, std::vector<double>, KO_Traits::StoringMatrix>; 

template <VALID_ERR_RET valid_err_ret>
using results_t = std::conditional<valid_err_ret,results_err_t,results_no_err_t>::type;



#endif /*KO_TRAITS_HPP*/