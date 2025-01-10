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

#ifndef STRAT_CVPPC_HPP
#define STRAT_CVPPC_HPP

#include <iostream>
#include <functional>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <concepts>
#include <memory>
#include <utility>
#include <type_traits>


#include "traits_ko.hpp"


/*!
* @file strategy_cv.hpp
* @brief Contains the class for creating the split training/validation set
* @author Andrea Enrico Franzoni
*/


/*!
* Doing tag dispatching for the correct way of splitting training/validation set.
* @tparam cv_strat: template parameter for the splitting strategy
*/
template <CV_STRAT cv_strat>
using CV_STRAT_T = std::integral_constant<CV_STRAT, cv_strat>;


/*!
* Type for splitting strategy (a vector that contains pairs: first element column indices of training set, second of validation set)
*/
using cv_strategy_t     = std::vector<std::pair<std::vector<int>,std::vector<int>>>;
/*!
* Type for single cv iter (a pair that identifies a specific training and valdiation set)
*/
using iter_cv_t         = std::pair<std::vector<int>,std::vector<int>>;
/*!
* Type for training and validation sets
*/
using train_valid_set_t = std::pair<KO_Traits::StoringMatrix,KO_Traits::StoringMatrix>;


/*!
* @class cv_strategy
* @brief Template class for creating training/validation split according to a specific strategy.
* @tparam cv_strat is the splitting training/validation strategy
*/
template<CV_STRAT cv_strat>
class cv_strategy
{
private:

  /*!
  * @brief Creating the training/validation split according to augmenting window strategy
  */
  void train_validation_set_strategy(int min_dim_ts, int max_dim_ts, CV_STRAT_T<CV_STRAT::AUGMENTING_WINDOW>);
  
  /*!Training/validation splitting*/
  cv_strategy_t m_strategy;
  
  /*!
  * @brief For a fixed given split training/validation according to augmenting window strategy, returns the two sets
  * @param data matrix containing the fts
  * @param strat a given split training/validation
  */
  train_valid_set_t train_validation_set(const KO_Traits::StoringMatrix &data, const iter_cv_t &strat, CV_STRAT_T<CV_STRAT::AUGMENTING_WINDOW>) const;
  
public:
  
  /*!
  * @brief Constructor taking minimum and maximum dimension of the training set
  * @param min_dim_ts minimum dimension of the training set
  * @param max_dim_ts maximum dimension of the training set
  */
  cv_strategy(int min_dim_ts, int max_dim_ts) { this->train_validation_set_strategy(min_dim_ts,max_dim_ts);}
  
  /*!
  * @brief Getter for the splitting training/validation
  * @return the private m_strategy
  */
  inline cv_strategy_t strategy() const {return m_strategy;}
  
  /*!
  * @brief Creating the training/validation split. Tag-dispacther.
  * @param min_dim_ts minimum dimension of the training set
  * @param max_dim_ts maximum dimension of the training set
  * @details the method is called in the constructor, and update the private member 'm_strategy'
  */
  void train_validation_set_strategy(int min_dim_ts, int max_dim_ts) { return train_validation_set_strategy(min_dim_ts, max_dim_ts, CV_STRAT_T<cv_strat>{});};

  /*!
  * @brief For a fixed given split training/validation according to augmenting window strategy, returns the two sets. Tag-dispacther.
  * @param data matrix containing the fts
  * @param strat a given split training/validation
  */
  train_valid_set_t train_validation_set(const KO_Traits::StoringMatrix &data, const iter_cv_t &strat) const { return train_validation_set(data, strat, CV_STRAT_T<cv_strat>{});};
  
};


#include "strategy_cv_imp.hpp"

#endif  //STRAT_CVPPC_HPP