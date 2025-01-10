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

#ifndef KO_UTILS_HPP
#define KO_UTILS_HPP

#include "traits_ko.hpp"
#include <limits>

/*!
* @file utils.hpp
* @brief Contains helpers functions
* @author Andrea Enrico Franzoni
*/


/*!
* @brief Function to dispacth the error type depending on if cv is performed on one or two parameters, using visitor design pattern
* @param var errors, stored in variant able to represent a vector or a matrix
* @return a list containing the error stored according to the dimension of the cv input space
*/
Rcpp::List 
valid_err_disp
(const valid_err_variant& var) 
{

  //visitor design pattern to dispatch the error correctly 
  return std::visit([](auto&& arg) -> Rcpp::List 
  {
    // retrieving the type of the variant object
    using T = std::decay_t<decltype(arg)>;
    
    // cv single param: errors in a vector -> only one vector back
    if constexpr (std::is_same_v<T, std::vector<double>>) 
    {
      return Rcpp::List::create(Rcpp::Named("Errors") = Rcpp::wrap(arg));  // NumericVector
    }
    // cv double param: errors in a vector of vector -> list of vector(one for every alpha) back
    else if constexpr (std::is_same_v<T, std::vector<std::vector<double>>>) 
    {
      Rcpp::List res(arg.size());
      for (size_t i = 0; i < arg.size(); ++i) 
      {
        res[i] = Rcpp::wrap(arg[i]);  
      }
      
      return Rcpp::List::create(Rcpp::Named("Errors") = res);              // List of NumericVector
    } 
    // check for the corret variant
    else 
    {
      Rcpp::stop("Wrong type in variant!");
    }
    }, var);
}


/*!
* @brief Function to map a matrix into a column vector, column by column
* @param mat matrix to be mapped
* @return a vector over which the matrix has been mapped
*/
KO_Traits::StoringVector
from_matrix_to_col(const KO_Traits::StoringMatrix &mat)
{
  return KO_Traits::StoringVector(mat.reshaped());
}


/*!
* @brief Function to map a column vector into a matrix, column by column
* @param col vector to be mapped
* @param rows number of rows of the matrix
* @param cols number of columns of the matrix
* @return a matrix over which the vector has been mapped
*/
KO_Traits::StoringMatrix
from_col_to_matrix(const KO_Traits::StoringVector &col, int rows, int cols)
{
  return Eigen::Map<const KO_Traits::StoringMatrix>(col.data(),rows,cols);
}


//function to add nans in the result
//pred is the vector with the prediction, row ret the vector containing which points actually have evals, complete size is the size of all the points
/*!
* @brief Function to add dummy NaNs to the curve/surface
* @param pred vector where NaNs have to be added
* @param row_ret vector containing the actual rows that are not dummy NaNs 
* @param complete_size size of the returning vector considering also the dummy NaNs
* @return a vector with dummy NaNs in the requested positions
*/
KO_Traits::StoringVector
add_nans_vec(const KO_Traits::StoringVector &pred, const std::vector<int> &row_ret, int complete_size)
{
  if(row_ret.size()==0){return pred;}
  KO_Traits::StoringVector pred_comp(complete_size);
  
  pred_comp.setConstant(std::numeric_limits<double>::quiet_NaN());
  
  // putting values where they actually are
  int counter=0;
  std::for_each(row_ret.cbegin(),row_ret.cend(),[&pred_comp,&pred,&counter](int el){pred_comp[el]=pred[counter]; counter++;});
  
  return pred_comp;
}

#endif  //KO_UTILS_HPP