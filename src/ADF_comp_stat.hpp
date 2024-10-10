#ifndef ADF_PPC_STAT_HPP
#define ADF_PPC_STAT_HPP

#include <array>

namespace ADF_util
{

//using a functor and the template parameter as policy
//computing the test statistic for the adf
template <class LAG_policy> 
class adf_stat
{
public:
  
  double
  operator()(const KO_Traits::StoringVector &x, const KO_Traits::StoringMatrix &z) const
  {
    //test statistic of ADF, computed according to the fact that there is lag order or not (policy via template)
    std::array<double,2> values = lag_order(x,z);
    return values[0]/values[1];
  }
  
private:
  LAG_policy lag_order;
};

}   //end namespace ADF_util

#endif  //ADF_PPC_STAT_HPP