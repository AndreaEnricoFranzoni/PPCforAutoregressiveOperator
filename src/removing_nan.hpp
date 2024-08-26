#ifndef KO_REMOVE_NAN_HPP
#define KO_REMOVE_NAN_HPP

#include <Eigen/Dense>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <functional>
#include <vector>
#include <type_traits>

#include "wrapper_params.hpp"
#include "default_parameters.hpp"
#include "KO_Traits.hpp"


//to do tag dispatching for the correct way of removing nans
template <DEF_PARAMS_PPC::MA_type MA_t>
using MAT = std::integral_constant<DEF_PARAMS_PPC::MA_type, MA_t>;


namespace REM_NAN_PPC           //utilities to actually remove nans
{

template<typename T,DEF_PARAMS_PPC::MA_type MA_t>
class removing_nan
{
public:
  
  //constructor: takes the matrix from which there is need to leave out the nans
  removing_nan(KO_Traits::StoringMatrix&& data)
    :
    m_data{std::forward<KO_Traits::StoringMatrix>(data)}
    {}
  
  //getter for m_data
  inline KO_Traits::StoringMatrix data() const {return m_data;};
  
  //how to leave out the nans in each row: tag dispatching is performed
  void row_removal(Eigen::Block<Eigen::Matrix<T,-1,-1>,1>& row) { return row_removal(row, MAT<MA_t>{});};
  
  //for each row: substitute the nans
  inline void remove_nan(){ for(auto row : m_data.rowwise()){row_removal(row);}};

  
private:
  KO_Traits::StoringMatrix m_data;
  std::size_t m_first_ti;
  std::size_t m_last_ti;
  std::vector<double> m_weights;
  
  void row_removal(Eigen::Block<Eigen::Matrix<T,-1,-1>,1>& row, MAT<DEF_PARAMS_PPC::MA_type::EMA>);
  void row_removal(Eigen::Block<Eigen::Matrix<T,-1,-1>,1>& row, MAT<DEF_PARAMS_PPC::MA_type::WMA>);
  void row_removal(Eigen::Block<Eigen::Matrix<T,-1,-1>,1>& row, MAT<DEF_PARAMS_PPC::MA_type::SMA>);
  void row_removal(Eigen::Block<Eigen::Matrix<T,-1,-1>,1>& row, MAT<DEF_PARAMS_PPC::MA_type::MR>);
  void row_removal(Eigen::Block<Eigen::Matrix<T,-1,-1>,1>& row, MAT<DEF_PARAMS_PPC::MA_type::ZR>);
  
};

} //end namespace REM_NAN_PPC

#include "removing_nan_imp.hpp"

#endif /*KO_REMOVE_NAN_HPP*/