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

#ifndef KO_READ_DATA_HPP
#define KO_READ_DATA_HPP

#include <RcppEigen.h>

#include <algorithm>
#include <iterator>
#include <numeric>
#include <cmath>
#include <utility>
#include <vector>
#include <iostream>

#include "traits_ko.hpp"
#include "removing_nan.hpp"
#include "parameters_wrapper.hpp"



/*!
* @file data_reader.hpp
* @brief Contains the method to read data
* @author Andrea Enrico Franzoni
*/


/*!
* @brief Function that reads data from R containers, substitute actual NaNs, handle dummy NaNs, removing them but saving their position, and wraps it into C++ objects
* @param X Rcpp::NumericMatrix as passed in input to the R-interfaced function
* @param MA_t how to handle not-dummy NaNs (if substituting with pointwise fts mean or 0)
* @return a pair containing the mapped matrix and a vector with the positions of the dummy NaNs (rows of the original matrix)
* @note Depends on RcppEigen for interfacing with R containers
*/
//
// [[Rcpp::depends(RcppEigen)]]
template<typename T> 
std::pair<KO_Traits::StoringMatrix,std::vector<int>>
reader_data(Rcpp::NumericMatrix X,
            REM_NAN MA_t)
{
  //taking the dimensions: n_row is the number of time series, n_col is the number of time istants
  int n_row = X.nrow();
  int n_col = X.ncol();
  
  //needed only if the method is used to translate R containers into a coherent input matrix for surface's PPCKO verson. NOT usable by the user
  if(MA_t == REM_NAN::NR)
  {
    KO_Traits::StoringMatrix x = Eigen::Map<KO_Traits::StoringMatrix>(X.begin(),n_row,n_col);
    std::vector<int> row_removed;
    return std::make_pair(x,row_removed);
  }
  

  //removing rows of all dummy NaNs
  std::set<int> rows_to_be_retained = rows_entire_NaNs(X);      //rows to be retained (position wrt the origianl matrix)
  
  std::vector<int> temp;
  temp.resize(X.nrow());
  std::iota(temp.begin(),temp.end(),static_cast<int>(0));     
  std::set<int> all_rows(temp.begin(),temp.end());
  temp.clear();
  
  std::set<int> rows_to_be_removed;                             //rows to be removed (position wrt the origianl matrix)
  std::set_difference(all_rows.begin(),all_rows.end(),
                      rows_to_be_retained.begin(),rows_to_be_retained.end(),
                      std::inserter(rows_to_be_removed,rows_to_be_removed.begin()));
  all_rows.clear();
  
  //remove rows of all dummy NaNs
  Rcpp::NumericMatrix X_clean = removing_NaNS_entire_rows(X,rows_to_be_retained);         //matrix without dummy NaNs
  //mapping the position of the rows to be retained into a std::vector from a std::set
  std::vector<int> rows_retained(rows_to_be_retained.begin(),rows_to_be_retained.end());  //AS IF INDECES START FROM 0
  
  
  //Eigen::Map to map data into KO_Traits::StoringMatrix (Eigen::MatrixXd)  (!!to be modified if the trait is not an Eigen object anymore!!)
  KO_Traits::StoringMatrix x = Eigen::Map<KO_Traits::StoringMatrix>(X_clean.begin(),rows_to_be_retained.size(),n_col);
  
  //check if there are other NaNs (NaNs due to missed measurements, not dummy)
  auto check_nan = std::find_if(x.reshaped().cbegin(),x.reshaped().cend(),[](T el){return std::isnan(el);});
  
  //if there are nans: remove them
  if (check_nan!=x.reshaped().end())
  {
    if(MA_t == REM_NAN::MR)     //replacing nans with the mean
    {
      removing_nan<T,REM_NAN::MR> data_clean(std::move(x));
      data_clean.remove_nan();
      return std::make_pair(data_clean.data(),rows_retained);       
    }
    if(MA_t == REM_NAN::ZR)     //replacing nans with 0s
    {
      removing_nan<T,REM_NAN::ZR> data_clean(std::move(x));
      data_clean.remove_nan();
      return std::make_pair(data_clean.data(),rows_retained);
    }
  }
  
  return std::make_pair(x,rows_retained);
}

#endif /*KO_READ_DATA_HPP*/