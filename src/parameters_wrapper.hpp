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


//utilities to wrap the input parameters of the R function



//to print which version is used
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


//check that threshold_ppc is in the correct range
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



//check that alpha is in the correct range
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



//check that k is in the correct range
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



//check that the solver is passed correctly
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


//to wrap the vector of alphas passed in input
inline
std::vector<double>
wrap_alpha_vec(Rcpp::Nullable<Rcpp::NumericVector> alpha_vec)
{
  //if no alpha is given, the default value for the alphas is given
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



//to wrap the vector of k passed in input
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
    
    //sorting into ascending order the alphas to be coherent during the algorithm
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



//to wrap the points in which the domain has been discretized
inline
std::vector<double>
wrap_disc_ev(Rcpp::Nullable<Rcpp::NumericVector> disc_ev, double a, double b, int dim)    //dim: row of x
{ 
  if(a>=b)
  {
    std::string error_message1 = "Left extreme of the domain has to be smaller than the right one";
    throw std::invalid_argument(error_message1);
  }
  
  if(disc_ev.isNull())
  {
    Geometry::Domain1D domain_func_data(a,b);
    Geometry::Mesh1D   grid_func_data(domain_func_data,dim-static_cast<int>(1));
    
    return grid_func_data.nodes();
  }
  
  std::vector<double> disc_ev_points = Rcpp::as<std::vector<double>>(disc_ev);
  std::sort(disc_ev_points.begin(),disc_ev_points.end());
  
  if(disc_ev_points[0] < a || disc_ev_points.back() > b)
  {
    std::string error_message2 = "The points in which there are the discrete evaluations of the functiona data have to in the domain (" + std::to_string(a) + "," + std::to_string(b) + ")";
    throw std::invalid_argument(error_message2);
  }
  
  if(disc_ev_points.size()!=dim)
  {
    std::string error_message3 = "In the grid are needed " + std::to_string(dim) + " points";
    throw std::invalid_argument(error_message3);
  }
  
  return disc_ev_points;
}

//to wrap the min a max dimension of ts
inline
std::pair<int,int>
wrap_sizes_set_CV(Rcpp::Nullable<int> min_size_ts, Rcpp::Nullable<int> max_size_ts, int number_time_instants)    //dim: row of x
{ 
  int min_dim_ts = min_size_ts.isNull() ? static_cast<int>(std::ceil(static_cast<double>(number_time_instants)/static_cast<double>(2)))  : Rcpp::as<int>(min_size_ts);
  int max_dim_ts = max_size_ts.isNull() ? number_time_instants  : (Rcpp::as<int>(max_size_ts)+1); 
  
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


//to wrap the number of threads
inline
int
wrap_num_thread(Rcpp::Nullable<int> num_threads)
{
#ifndef _OPENMP
  return 1;
#else
  
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
enum REM_NAN
{ 
  NR = 0,      //not replacing NaN
  MR = 1,      //replacing nans with mean (easily changes the mean of the distribution)
  ZR = 2,      //replacing nans with 0s (easily changes the sd of the distribution)
};


//reads the input string an gives back the correct value of the enumerator for replacing nans
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
