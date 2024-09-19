#include "removing_nan.hpp"

//ExponentailMovingAverage
template<typename T,DEF_PARAMS_PPC::MA_type MA_t>
void
REM_NAN_PPC::removing_nan<T,MA_t>::row_removal(Eigen::Block<Eigen::Matrix<T,-1,-1>,1>& row, MAT<DEF_PARAMS_PPC::MA_type::EMA>)
{
  std::cout << "I am using EMA with td" << std::endl;
  std::replace_if(row.begin(),row.end(),[](T el){return isnan(el);},static_cast<T>(0));
}  



//WeightedMovingAverage
template<typename T,DEF_PARAMS_PPC::MA_type MA_t>
void
REM_NAN_PPC::removing_nan<T,MA_t>::row_removal(Eigen::Block<Eigen::Matrix<T,-1,-1>,1>& row, MAT<DEF_PARAMS_PPC::MA_type::WMA>)
{
  std::cout << "I am using WMA with td" << std::endl;
  std::replace_if(row.begin(),row.end(),[](T el){return isnan(el);},static_cast<T>(0));
}  



//SimpleMovingAverage
template<typename T,DEF_PARAMS_PPC::MA_type MA_t>
void
REM_NAN_PPC::removing_nan<T,MA_t>::row_removal(Eigen::Block<Eigen::Matrix<T,-1,-1>,1>& row, MAT<DEF_PARAMS_PPC::MA_type::SMA>)
{
  std::cout << "I am using SMA with td" << std::endl;
  std::replace_if(row.begin(),row.end(),[](T el){return isnan(el);},static_cast<T>(0));
}  



//MeanReplace
template<typename T,DEF_PARAMS_PPC::MA_type MA_t>
void
REM_NAN_PPC::removing_nan<T,MA_t>::row_removal(Eigen::Block<Eigen::Matrix<T,-1,-1>,1>& row, MAT<DEF_PARAMS_PPC::MA_type::MR>)
{
  std::cout << "I am using MR with td" << std::endl;
  const auto n_not_nan = std::count_if(row.begin(),row.end(),[](T el){return !isnan(el);});
  std::vector<T> el_not_nan;
  el_not_nan.reserve(n_not_nan);
  std::copy_if(row.begin(),row.end(),std::back_inserter(el_not_nan),[](T el){return !isnan(el);});
  const T mean = std::accumulate(el_not_nan.begin(),el_not_nan.end(),static_cast<T>(0),std::plus{})/static_cast<T>(n_not_nan);
  el_not_nan.clear();
  std::replace_if(row.begin(),row.end(),[](T el){return isnan(el);},mean);
}  



//ZerosReplace
template<typename T,DEF_PARAMS_PPC::MA_type MA_t>
void
REM_NAN_PPC::removing_nan<T,MA_t>::row_removal(Eigen::Block<Eigen::Matrix<T,-1,-1>,1>& row, MAT<DEF_PARAMS_PPC::MA_type::ZR>)
{
  std::cout << "I am using ZR with td" << std::endl;
  std::replace_if(row.begin(),row.end(),[](T el){return isnan(el);},static_cast<T>(0)); 
}   