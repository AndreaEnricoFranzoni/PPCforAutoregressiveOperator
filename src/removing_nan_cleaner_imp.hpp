#include "removing_nan.hpp"

template<typename T,DEF_PARAMS_PPC::MA_type MA_t>
KO_Traits::StoringMatrix
REM_NAN_PPC::removing_nan<T,MA_t>::final_clean()
const
{ 
  std::vector<int> temp;
  temp.resize(m_m);
  std::iota(temp.begin(),temp.end(),static_cast<int>(1));     
  
  //all the rows
  std::set<int> all_row(temp.begin(),temp.end());
  temp.clear();
  
  //how many rows have to be removed
  auto n_row_rem = m_row_to_be_removed.size();
  
  //if there are only NaNs: error
  if(n_row_rem==m_m)
  {
    std::string error_message = "Input data are all NaNs";
    throw std::invalid_argument(error_message);  
  }

  //rows retained: only the ones that are not removed
  std::set<int> row_ret;
  std::set_difference(all_row.begin(),all_row.end(),
                      m_row_to_be_removed.begin(),m_row_to_be_removed.end(),
                      std::inserter(row_ret,row_ret.begin()));
  
  //remove the rows of all NaNs
  KO_Traits::StoringMatrix data_cleaned(m_m - n_row_rem, m_n);
  int counter_row_ret = 0;
  std::for_each(row_ret.begin(),row_ret.end(),[this,&data_cleaned,&counter_row_ret](int el){data_cleaned.row(counter_row_ret)=this->m_data.row(el-1); counter_row_ret++;});

  return data_cleaned;
}