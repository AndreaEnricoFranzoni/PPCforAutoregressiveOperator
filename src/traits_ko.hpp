#ifndef KO_TRAITS_HPP
#define KO_TRAITS_HPP

#include <Eigen/Dense>
#include <tuple>
#include <vector>
#include <array>
#include <variant>
#include <type_traits>
#include <cmath>
#include <string>


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


//algo implemented
struct CV_algo
{
  static constexpr std::string CV1 = "NoCV";
  static constexpr std::string CV2 = "CV_alpha";
  static constexpr std::string CV3 = "CV_k";
  static constexpr std::string CV4 = "CV";
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

//results from the algo
//        EigenVector with one-step ahead prediction:   Y
//        regularization parameter used                 Y
//        number of PPCs                                Y
//        scores along the PPCs                         Y
//        explanatory power of the PPCs                 Y
//        a_i (directions)                              Y
//        b_i (weights)                                 Y
//        id_CV (which type of algorithm has been used)
//        left extreme of the domain
//        right extreme of the domain
//        grid for the discrete evaluations
//        errors (eventually)                           Y
using results_err_t = std::tuple<KO_Traits::StoringVector, double, int, std::vector<double>, std::vector<double>, KO_Traits::StoringMatrix, KO_Traits::StoringMatrix, valid_err_variant>; 
using results_no_err_t = std::tuple<KO_Traits::StoringVector, double, int, std::vector<double>, std::vector<double>, KO_Traits::StoringMatrix, KO_Traits::StoringMatrix>; 

// if errors are saved or not
template <VALID_ERR_RET valid_err_ret>
using results_t = typename std::conditional<valid_err_ret,results_err_t,results_no_err_t>::type;



#endif /*KO_TRAITS_HPP*/
