#include "removing_nan.hpp"

//MeanReplace
template<typename T,REM_NAN MA_t>
void
removing_nan<T,MA_t>::row_removal(Eigen::Block<Eigen::Matrix<T,-1,-1>,1>& row, MAT<REM_NAN::MR>)
{
  const auto n_not_nan = std::count_if(row.begin(),row.end(),[](T el){return !isnan(el);});
  std::vector<T> el_not_nan;
  el_not_nan.reserve(n_not_nan);
  std::copy_if(row.begin(),row.end(),std::back_inserter(el_not_nan),[](T el){return !isnan(el);});
  const T mean = std::accumulate(el_not_nan.begin(),el_not_nan.end(),static_cast<T>(0),std::plus{})/static_cast<T>(n_not_nan);
  el_not_nan.clear();
  std::replace_if(row.begin(),row.end(),[](T el){return isnan(el);},mean);
}  



//ZerosReplace
template<typename T,REM_NAN MA_t>
void
removing_nan<T,MA_t>::row_removal(Eigen::Block<Eigen::Matrix<T,-1,-1>,1>& row, MAT<REM_NAN::ZR>)
{
  std::replace_if(row.begin(),row.end(),[](T el){return isnan(el);},static_cast<T>(0)); 
} 