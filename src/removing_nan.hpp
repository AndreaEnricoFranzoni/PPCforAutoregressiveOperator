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
    { 
      m_m = m_data.rows();
      m_n = m_data.cols();
    }
  
  //getter for m_data
  inline KO_Traits::StoringMatrix data() const {return m_data;};
  
  //removing rows with too many nans
  KO_Traits::StoringMatrix final_clean() const;
  
  //how many nans to actually remove an entire row?
  inline bool
  const
  cond_rem_row(Eigen::Block<Eigen::Matrix<T,-1,-1>,1>& row)
  { 
    //condition for removing a row since there are too many nans: TODO: change this?
    //all nans in the row
    return std::count_if(row.begin(),row.end(),[](T el){return isnan(el);}) == m_n  ?  true  :  false;
  };
  
  //how to leave out the nans in each row: tag dispatching is performed to substitute them correctly
  void row_removal(Eigen::Block<Eigen::Matrix<T,-1,-1>,1>& row) { return row_removal(row, MAT<MA_t>{});};
  
  //for each row: if not too many nans: substitute them; else, the row has to be removed
  inline void remove_nan()
  { int counter_row = 1;
    for(auto row : m_data.rowwise()){ if(cond_rem_row(row)){m_row_to_be_removed.insert(counter_row);} else{row_removal(row);};   counter_row++;} 
    if(!m_row_to_be_removed.empty()){ m_data = final_clean();}   
  };
  
  
private:
  KO_Traits::StoringMatrix m_data;
  std::size_t m_m;
  std::size_t m_n;
  std::set<int> m_row_to_be_removed;
  std::size_t m_first_ti;
  std::size_t m_last_ti;
  std::vector<double> m_weights;
  

  void row_removal(Eigen::Block<Eigen::Matrix<T,-1,-1>,1>& row, MAT<DEF_PARAMS_PPC::MA_type::MR>);
  void row_removal(Eigen::Block<Eigen::Matrix<T,-1,-1>,1>& row, MAT<DEF_PARAMS_PPC::MA_type::ZR>);
  
};

} //end namespace REM_NAN_PPC

#include "removing_nan_imp.hpp"
#include "removing_nan_cleaner_imp.hpp"

#endif /*KO_REMOVE_NAN_HPP*/