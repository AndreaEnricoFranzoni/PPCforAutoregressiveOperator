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

/*!
* @file traits_ko.hpp
* @brief Contains customized types and enumerator for customized template parameters, exploited in the algorithm
* @author Andrea Enrico Franzoni
*/



/*!
* @struct KO_Traits
* @brief Contains the customized types for fts, covariances, PPCs, etc...
* @details Data are stored in dynamic matrices (easily very big dimensions) of doubles
*/
struct KO_Traits
{
public:
  
  using StoringMatrix = Eigen::MatrixXd;  ///< Matrix data structure.
  
  using StoringVector = Eigen::VectorXd;  ///< Vector data structure.
  
  using StoringArray  = Eigen::ArrayXd;   ///< Array data structure: more efficient for coefficient-wise operations.

};


/*!
* @struct CV_algo
* @brief Contains PPCKO versions implemented
*/
struct CV_algo
{
  static constexpr std::string CV1 = "NoCV";        ///< No cv for parameters.
  static constexpr std::string CV2 = "CV_alpha";    ///< Cv for regularization parameter.
  static constexpr std::string CV3 = "CV_k";        ///< Cv for number of retained PPCs.
  static constexpr std::string CV4 = "CV";          ///< Cv for both regularization parameter and number of retained PPCs.
};


/*!
* @enum SOLVER
* @brief The available solvers for PPCKO algorithm
*/
enum SOLVER
{
  ex_solver  = 0,      ///< Inverted square root regularzied covariance and retrieving PPCs from phi
  gep_solver = 1,      ///< Using GEP to avoid to avoid inverted square root
};


/*!
* @enum K_IMP
* @brief If the number of PPCs k is a parameter 
*/
enum K_IMP
{
  NO  = 0,     ///< k is not passed as parameter, has to be found using cumulative explanatory power
  YES = 1,     ///< k is already known (it will be fixed if doing cv on it) 
};


/*!
* @enum VALID_ERR_RET
* @brief If validation error has to be stored and returned 
*/
enum VALID_ERR_RET
{
  NO_err   = 0,     ///< Validation errors are not stored and not returned (memory saving)
  YES_err  = 1,     ///< Validation errors are stored and returned
};


/*!
* @enum CV_STRAT
* @brief Strategy for training/validation splitting during cv
*/
enum CV_STRAT
{
  AUGMENTING_WINDOW = 0,  ///< Fixing an instant: training set are all the instant up to it, validation the next one. At each iteration: the previous validation set is inglobated into the previous training set for the new one. Validation set is shifted by one
};


/*!
* @enum CV_ERR_EVAL
* @brief How to compute validation errors during cv
*/
enum CV_ERR_EVAL
{
  MSE = 0,  ///< Estimate of the L2 norm loss
};


/*!
* Types for the errors: variant is used (for cv on both parameter a matrix is returned, a vector otherwise)
*/
using valid_err_cv_1_t = std::vector<double>;
using valid_err_cv_2_t = std::vector<std::vector<double>>;
using valid_err_variant = std::variant<valid_err_cv_1_t,valid_err_cv_2_t>;


/*!
* Types for the results: tuple is exploited, dimension and types depends on if valdiation errors are returned and which ones eventually
*/
using results_err_t = std::tuple<KO_Traits::StoringVector, double, int, std::vector<double>, std::vector<double>, KO_Traits::StoringMatrix, KO_Traits::StoringMatrix, std::vector<std::array<double,2>>, KO_Traits::StoringArray, valid_err_variant>; 
using results_no_err_t = std::tuple<KO_Traits::StoringVector, double, int, std::vector<double>, std::vector<double>, KO_Traits::StoringMatrix, KO_Traits::StoringMatrix, std::vector<std::array<double,2>>, KO_Traits::StoringArray>; 

/*!
* Type for the returning error: depending on
* @param valid_err_ret: enumerator VALID_ERR_RET: if errors have to be returned on not
*/
template <VALID_ERR_RET valid_err_ret>
using results_t = typename std::conditional<valid_err_ret,results_err_t,results_no_err_t>::type;  //if errors are returned: the tuple will have one more element



#endif /*KO_TRAITS_HPP*/