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

#include "CV.hpp"


/*!
* @file CV_err_valid_eval.hpp
* @brief Implementation of tag-dispatched function for evaluating the validation error
* @author Andrea Enrico Franzoni
*/


/*!
* @brief Error evaluation between a prediction on validation set and validation set (in a fixed cv iteration). L2 norm estimate of the error.
* @param pred prediction on validation set
* @param valid validation set
* @param number_threads number of threads for OMP
* @details 'CV_ERR_EVAL::MSE' dispatch.
*/
template< class D, CV_STRAT cv_strat, CV_ERR_EVAL err_eval, K_IMP k_imp, VALID_ERR_RET valid_err_ret >
double
CV_base<D,cv_strat,err_eval,k_imp,valid_err_ret>::err_valid_set_eval(const KO_Traits::StoringVector &pred, const KO_Traits::StoringVector &valid, int number_threads, ERR_EVAL_T<CV_ERR_EVAL::MSE>)
const
{
  //using mse between predicted and validation
  return mse<double>( pred.array() - valid.array(), number_threads );
}