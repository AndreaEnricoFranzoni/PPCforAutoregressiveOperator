// Copyright (c) 2024 Andrea Enrico Franzoni (andreaenrico.franzoni@gmail.com)
//
// This file is part of PPCKO
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of PPCKO and associated documentation files (the PPCKO software), to deal
// PPCKO without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of PPCKO, and to permit persons to whom PPCKO is
// furnished to do so, subject to the following conditions:
//
// PPCKO IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH PPCKO OR THE USE OR OTHER DEALINGS IN
// PPCKO.

#ifndef KO_WRAP_PARAMS_HPP
#define KO_WRAP_PARAMS_HPP

#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <utility>
#include <string>
#include <stdexcept>
#include "traits_ko.hpp"

#include "mesh.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif


/*!
* @file parameters_wrapper.hpp
* @brief Contains methods to check and wrap R-inputs into PPCKO-coherent ones.
* @author Andrea Enrico Franzoni
*/



/*!
* @brief Creating a string to print on the screen which PPCKO version is being used
* @param id_cv string indicating the version used
* @return a string to be printed
*/
inline 
std::string
wrap_string_CV_to_be_printed(const std::string &id_cv)
{ 
  if(id_cv==CV_algo::CV1){  return "no cross-validation";}
  if(id_cv==CV_algo::CV2){  return "cross validation on regularization parameter";}
  if(id_cv==CV_algo::CV3){  return "cross validation on number of PPCs";}
  if(id_cv==CV_algo::CV4){  return "cross validation on both regularization parameter and number of PPCs";}
  else
  {
    std::string error_message = "Wrong input string";
    throw std::invalid_argument(error_message);
  }
}



/*!
* @brief Check if 'threshold_ppc' input is between 0 and 1. Eventually, raises and error.
* @param threshold_ppc requested explanatory power for the predictor
*/
inline
void
check_threshold_ppc(const double &threshold_ppc)
{
  if(threshold_ppc<=0 || threshold_ppc>=1)
  {
    std::string error_message = "threshold_ppc has to be in (0,1)";
    throw std::invalid_argument(error_message);
  }
}



/*!
* @brief Check if 'alpha' input is greater than 0. Eventually, raises and error.
* @param alpha regularization parameter
*/
inline
void
check_alpha(const double &alpha)
{
  if( alpha<= 0 )
  {
    std::string error_message = "alpha has to be a positive real number";
    throw std::invalid_argument(error_message);
  }
}



/*!
* @brief Check if 'k' input is an integer between 0 and the number of available evaluations of the functional object. Eventually, raises and error.
* @param k number of PPCs passed as parameter
* @param max_k the maximum value for 'k' (number of available evaluations of the functional object)
*/
inline
void
check_k(const int &k, const int &max_k)
{
  if( k < 0 )
  {
    std::string error_message = "k has to be a positive integer or 0";
    throw std::invalid_argument(error_message);
  }
  if( k > max_k )
  {
    std::string error_message = "k has to be lower than the maximum number of PPCs";
    throw std::invalid_argument(error_message);
  }
}



/*!
* @brief Check if, if using 'gep_solver', the number of PPCs is not retrieved through explanatory power criterion. Eventually, raises and error.
* @param solver_ex 'true' if using ex_solver
* @param id_cv which PPCKO version is used
* @param k the input parameter k
*/
inline
void
check_solver(bool solver_ex, const std::string &id_cv, int k)
{
  if(solver_ex == false)
  {
    if ( (k==0 && (id_cv==CV_algo::CV1 || id_cv==CV_algo::CV2))  )
    {
      std::string error_message = "GEP solver can be used only if the number of PPCs is imposed or retrieved through CV";
      throw std::invalid_argument(error_message);
    }
  }
}



/*!
* @brief Wrapping the R-vector representing the regularization parameter input space into a coherent C++ object, checking parameters consistency, eventually throwing an error, eventually sorting them in increasing order.
* @param alpha_vec Rcpp::Nullable<Rcpp::NumericVector> 
* @return the R-vector is shifted to a std::vector<double>
* @note If the input is 'NULL': a vector in logarithmic scale, from 10 to the power of -10 up to the power of 11 is generated
*/
inline
std::vector<double>
wrap_alpha_vec(Rcpp::Nullable<Rcpp::NumericVector> alpha_vec)
{
  //if no alphas are given, the default value for the alphas is generated
  if(alpha_vec.isNull())
  {
    std::vector<double> alphas;
    alphas.resize(21);
    
    std::iota(alphas.begin(),alphas.end(),static_cast<double>(-10));
    std::transform(alphas.begin(),alphas.end(),alphas.begin(),[](double el){return(pow(static_cast<double>(10),el));});
    
    return alphas;
  }
  
  std::vector<double> alphas = Rcpp::as<std::vector<double>>(alpha_vec);
  
  //sorting into ascending order the alphas to be coherent during the algorithm
  std::sort(alphas.begin(), alphas.end());
  
  if(alphas[0] <= 0)
  {
    std::string error_message = "Every alpha has to be a positive real number";
    throw std::invalid_argument(error_message);
  }
  
  return alphas;
}



/*!
* @brief Wrapping the R-vector representing the number of PPCs input space into a coherent C++ object, checking parameters consistency, eventually throwing an error, eventually sorting them in increasing order.
* @param k_vec Rcpp::Nullable<Rcpp::IntegerVector> 
* @param k_max the maximum number possible of retained PPCs
* @return the R-vector is shifted to a std::vector<int>
* @note If the input is 'NULL': a vector from 1 to 'k_max' is generated
*/
inline
std::vector<int>
wrap_k_vec(Rcpp::Nullable<Rcpp::IntegerVector> k_vec, int k_max)
{
  //if no k is given, the default is looking for all the possible directions (number of discrete evaluations innthe domain of the functional object)
    if(k_vec.isNull())
    {
      std::vector<int> k_s;
      k_s.resize(k_max);
      
      std::iota(k_s.begin(),k_s.end(),static_cast<int>(1));
    
      return k_s;
    }
    
    std::vector<int> k_s = Rcpp::as<std::vector<int>>(k_vec);
    
    //sorting into ascending order the ks to be coherent during the algorithm
    std::sort(k_s.begin(), k_s.end());
    
    //checking
    if(k_s[0] < 1)
    {
      std::string error_message1 = "k has to be at least 1";
      throw std::invalid_argument(error_message1);
    }
    if(k_s.back() > k_max)
    {
      std::string error_message2 = "k cannot be greater than the number of discrete evaluation of the functional object in the domain (" + std::to_string(k_max) + ")";
      throw std::invalid_argument(error_message2);
    }
    
    return k_s;
}



/*!
* @brief Wrapping the points over which the discrete evaluations of the functional object are available. Check consistency of domain extremes and passed points, eventualy throwing an error.
* @param disc_ev Rcpp::Nullable<Rcpp::NumericVector>  containing the domain points
* @param a left domain extreme
* @param b right domain extreme
* @param dim the number of discrete evaluations
* @return an std::vector<double>
* @note If the 'disc_ev' is 'NULL': an equally spaced grid from 'a' to 'b' is generated. For surfaces, the two dimensions are wrapped separately
*/
inline
std::vector<double>
wrap_disc_ev(Rcpp::Nullable<Rcpp::NumericVector> disc_ev, double a, double b, int dim)    //dim: row of x
{ 
  //check that domain extremes are consistent
  if(a>=b)
  {
    std::string error_message1 = "Left extreme of the domain has to be smaller than the right one";
    throw std::invalid_argument(error_message1);
  }
  
  //eventual default generation of the grid
  if(disc_ev.isNull())
  {
    Geometry::Domain1D domain_func_data(a,b);
    Geometry::Mesh1D   grid_func_data(domain_func_data,dim-static_cast<int>(1));
    
    return grid_func_data.nodes();
  }
  
  //sorting the abscissas values
  std::vector<double> disc_ev_points = Rcpp::as<std::vector<double>>(disc_ev);
  std::sort(disc_ev_points.begin(),disc_ev_points.end());
  
  //checking that the passed points are inside the domain
  if(disc_ev_points[0] < a || disc_ev_points.back() > b)
  {
    std::string error_message2 = "The points in which there are the discrete evaluations of the functiona data have to in the domain (" + std::to_string(a) + "," + std::to_string(b) + ")";
    throw std::invalid_argument(error_message2);
  }
  
  //checking that the dimension is ok (important when wrapping both the grid for surface case)
  if(disc_ev_points.size()!=dim)
  {
    std::string error_message3 = "In the grid are needed " + std::to_string(dim) + " points";
    throw std::invalid_argument(error_message3);
  }
  
  return disc_ev_points;
}



/*!
* @brief Wrapping the minimum and maximum dimension of the training set, checking their consitency, eventually raising an error.
* @param min_size_ts minimum dimension of the training set
* @param max_size_ts maximum dimension of the training set
* @param number_time_instants number of total time instants available
* @return a pair containig minimum and maximum dimension training set
* @note 
*   - 'min_dim_ts' has to be at least 2, but less than 'max_size_ts'. Default value: ceil of 'number_time_instants'/2
*   - 'max_size_ts' has to be greater than 'min_dim_ts' but less than number_time_instants. Default value: 'number_time_instants' (to be consistent with Eigen, look at the reference to better understand this choice)
* @see train_validation_set_strategy()
*/
inline
std::pair<int,int>
wrap_sizes_set_CV(Rcpp::Nullable<int> min_size_ts, Rcpp::Nullable<int> max_size_ts, int number_time_instants)    
{ 
  //eventual default value for min_size_ts: ceil number_time_instants/2
  int min_dim_ts = min_size_ts.isNull() ? static_cast<int>(std::ceil(static_cast<double>(number_time_instants)/static_cast<double>(2)))  : Rcpp::as<int>(min_size_ts);
  //eventual default value for max_size_ts: n
  int max_dim_ts = max_size_ts.isNull() ? number_time_instants  : (Rcpp::as<int>(max_size_ts)+1); 
  
  //checking that 
  if (!(min_dim_ts>1))
  {
    std::string error_message1 = "Min size of train set has to be at least 2";
    throw std::invalid_argument(error_message1);
  }

  if (!(max_dim_ts<=number_time_instants))
  {
    std::string error_message2 = "Max size of train set has to be at most the total number of time instants ("  + std::to_string(number_time_instants) + ") minus 1 (to leave room for the validation set)";
    throw std::invalid_argument(error_message2);
  }
  
  if(min_dim_ts >= max_dim_ts)
  {
    std::string error_message = "Min size of train set (" + std::to_string(min_dim_ts) + " has to be less than the max one ("  + std::to_string(max_dim_ts) + ")";
    throw std::invalid_argument(error_message);
  }

  return std::make_pair(min_dim_ts,max_dim_ts);
}



/*!
* @brief Wrapping the number of threads for OMP
* @param num_threads indicates how many threads to be used by multi-threading directives.
* @return the number of threads
* @details if omp is not included: will return 1. If not, a number going from 1 up to the maximum cores available by the machine used (default, or if the input is smaller than 1 or bigger than the maximum number of available cores)
* @note omp requested
*/
inline
int
wrap_num_thread(Rcpp::Nullable<int> num_threads)
{
#ifndef _OPENMP
  return 1;
#else
  
  //getting maximum number of cores in the machine
  int max_n_t = omp_get_num_procs();
  
  if(num_threads.isNull())
  {
    return max_n_t;
  }
  else
  {
    int n_t = Rcpp::as<int>(num_threads);
    if(n_t < 1 || n_t > max_n_t){  return max_n_t;}
    
    return n_t;
  }
#endif
}



//removing NaNs
/*!
* @enum REM_NAN
* @brief The available strategy for removing non-dummy NaNs
*/
enum REM_NAN
{ 
  NR = 0,      ///<  Not replacing NaN: not to be used by the user, necessary for handling dummy NaNs
  MR = 1,      ///< Replacing nans with mean (could change the mean of the distribution)
  ZR = 2,      ///< Replacing nans with 0s (could change the sd of the distribution)
};



/*!
* @brief Wrapping the strategy for handling non-dummy NaNs
* @param id_rem_nan string indicating the straegy for removing non-dummy NaNs
* @return the correpsonding value of 'REM_NAN' (default: 'MR')
*/
inline
REM_NAN
wrap_id_rem_nans(Rcpp::Nullable<std::string> id_rem_nan)
{
  if(id_rem_nan.isNull())
  { 
    return REM_NAN::MR;
  }
  if(Rcpp::as< std::string >(id_rem_nan) == "NO")
  {
    return REM_NAN::NR;
  }
  if(Rcpp::as< std::string >(id_rem_nan) == "MR")
  {
    return REM_NAN::MR;
  }
  if(Rcpp::as< std::string >(id_rem_nan) == "ZR")
  {
    return REM_NAN::ZR;
  }
  else
  {
    std::string error_message = "Wrong input string for handling NANs";
    throw std::invalid_argument(error_message);
  }
  
};


#endif  /*KO_WRAP_PARAMS_HPP*/