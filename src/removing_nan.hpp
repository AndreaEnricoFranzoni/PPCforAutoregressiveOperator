#ifndef KO_REMOVE_NAN_HPP
#define KO_REMOVE_NAN_HPP

#include <Eigen/Dense>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <functional>
#include <vector>
#include <set>
#include <type_traits>
#include <string>
#include <stdexcept>

#include "parameters_wrapper.hpp"
#include "traits_ko.hpp"



//to do tag dispatching for the correct way of removing nans
template <REM_NAN MA_t>
using MAT = std::integral_constant<REM_NAN, MA_t>;



template<typename T,REM_NAN MA_t>
class removing_nan
{
  
private:
  KO_Traits::StoringMatrix m_data;
  std::size_t m_m;
  std::size_t m_n;

  void row_removal(Eigen::Block<Eigen::Matrix<T,-1,-1>,1>& row, MAT<REM_NAN::MR>);
  void row_removal(Eigen::Block<Eigen::Matrix<T,-1,-1>,1>& row, MAT<REM_NAN::ZR>);
  
public:
  
  //constructor: takes the matrix from which there is need to leave out the nans
  template<typename STOR_OBJ>
  removing_nan(STOR_OBJ&& data)
    :
    m_data{std::forward<STOR_OBJ>(data)}
    { 
      m_m = m_data.rows();
      m_n = m_data.cols();
    }
  
  //getter for m_data
  inline KO_Traits::StoringMatrix data() const {return m_data;};
  
  //how to leave out the nans in each row: tag dispatching is performed to substitute them correctly
  void row_removal(Eigen::Block<Eigen::Matrix<T,-1,-1>,1>& row) { return row_removal(row, MAT<MA_t>{});};
  
  //for each row: if not too many nans: substitute them; else, the row has to be removed
  inline void remove_nan(){   for(auto row : m_data.rowwise()){   row_removal(row);} };
};


#include "removing_nan_imp.hpp"
#include "removing_nan_cleaner_imp.hpp"

#endif /*KO_REMOVE_NAN_HPP*/