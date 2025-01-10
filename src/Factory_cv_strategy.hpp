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

#ifndef CV_STRAT_FACTORY_HPP
#define CV_STRAT_FACTORY_HPP


#include <memory>
#include <utility>

#include "traits_ko.hpp"
#include "strategy_cv.hpp"



/*!
* @file Factory_cv_strategy.hpp
* @brief Factory for generating the training-validation splitting according to the strategy given
* @author Andrea Enrico Franzoni
*/



/*!
* @class Factory_cv_strat
* @brief Generating the training-validation splitting according to a given strategy
* @tparam cv_strat strategy for training/validation splitting
*/
template<CV_STRAT cv_strat> 
class Factory_cv_strat
{	
public:
  /*!
  * @brief Generating the training/validation splitting according to a given strategy
  * @tparam Args variadic template
  * @param args variadic number of arguments
  * @return a unique pointer to a cv_strategy<cv_strat> object
  * @details the method is static. std::unique_ptr provides a safer architecture for factory purposes
  */
  template<typename... Args>
  static 
  std::unique_ptr<cv_strategy<cv_strat>> 
  cv_strat_obj(Args&&...args)
    {
      //augmenting window
      if constexpr(cv_strat == CV_STRAT::AUGMENTING_WINDOW)
        { return std::make_unique<cv_strategy<cv_strat>>(std::forward<Args>(args)...);}
    }
};

#endif //CV_STRAT_FACTORY_HPP