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

#include "strategy_cv.hpp"


/*!
* @file strategy_cv_imp.hpp
* @brief Definition of methods for creating training/validation splitting and recovering corresponding training/validation sets
* @author Andrea Enrico Franzoni
*/


/*!
* @brief Creation of training/validation split accordinf to augmenting window strategy.
* @param min_dim_ts minimum dimension of training set
* @param max_dim_ts maximum dimension of training set
* @details 'CV_STRAT::AUGMENTING_WINDOW' dispatch. Modifying 'm_strategy' class private member. (Used as example to better understand @ref wrap_sizes_set_CV() "parameters_wrapper.hpp")
*/
template<CV_STRAT cv_strat>
void
cv_strategy<cv_strat>::train_validation_set_strategy(int min_dim_ts, int max_dim_ts, CV_STRAT_T<CV_STRAT::AUGMENTING_WINDOW>)
{
  std::size_t min_size_ts = static_cast<std::size_t>(min_dim_ts);
  std::size_t max_size_ts = static_cast<std::size_t>(max_dim_ts);
  
  //training set is defined as the instant from the beginning up to a given one
  //validation set is defined as the next one after the end of training set
  m_strategy.reserve(max_size_ts - min_size_ts);

  for(std::size_t i = min_size_ts; i < max_size_ts; ++i)
  { 
    //train set: i: Eigen takes the first i columns starting from 1
    std::vector<int> train_set;
    train_set.reserve(1);
    train_set.emplace_back(i);
    
    //validation set: i: Eigen takes the i-th column starting from 0
    std::vector<int> validation_set;
    validation_set.reserve(1);
    validation_set.emplace_back(i);
    
    m_strategy.emplace_back(std::make_pair(train_set,validation_set));
  }
}



/*!
* @brief Retaining a specific pair training and validation set given them as input.
* @param data matrix containing the fts
* @param strat a given pair training/validation set
* @return a pair: first element is the training set. Second element is the validation set.
* @details 'AUGMENTING_WINDOW' dispatch. Modifying 'm_strategy' class private member
*/
template<CV_STRAT cv_strat>
train_valid_set_t
cv_strategy<cv_strat>::train_validation_set(const KO_Traits::StoringMatrix &data, const iter_cv_t &strat, CV_STRAT_T<CV_STRAT::AUGMENTING_WINDOW>)
const
{
  return std::make_pair( data.leftCols(strat.first.front()), data.col(strat.second.front()) );
}