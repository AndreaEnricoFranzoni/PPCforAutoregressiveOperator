#ifndef ADF_PPC_HPP
#define ADF_PPC_HPP


#include <vector>
#include <algorithm>

#include "KO_Traits.hpp"
#include "ADF_comp_stat.hpp"
#include "ADF_comp_pvalue_util.hpp"



namespace PPC_util
{


//class to make ADF on a dataset such that each one of its row is a time serie
template<class LAG_policy>
class adf
{


private:
  KO_Traits::StoringMatrix m_x;               //data: each row (m) is a temporal series, and will be tested separately
  int m_k;                                    //if the lag order is used
  int m_k_used;
  int m_tot_time_instants;                    //number of time instants of the temporal serie
  std::vector<double> m_p_values;             //vector containing the pvalues of the m tests
 
 
public:
  adf(KO_Traits::StoringMatrix&& x, int k)      //constructor
   : m_x{std::forward<KO_Traits::StoringMatrix>(x)}, m_k(k), m_k_used(k+static_cast<int>(1)) 
   { m_tot_time_instants = m_x.cols(); } 
  
  /*!
   * Getter for m_p_values
   */
  inline std::vector<double> p_values() const {return m_p_values;};
  
  //function to evaluat one time step differences
  std::vector<double>      one_step_diff(const KO_Traits::StoringVector &ts)        const;
  KO_Traits::StoringMatrix embed(const KO_Traits::StoringVector &ts, int dimension) const;
  double                   statistic_eval(const KO_Traits::StoringVector &ts)       const;      //evaluation of the test statistic for the test on a single ts: it depends on k
  double                   p_value_eval(const KO_Traits::StoringVector &ts, const std::vector<double> &tableipl, int i)  const;                    //p_value evalution for the test on a single ts
  void                     test();
  
};
  
}   //end namespace PPC_util

#include "ADF_test_imp.hpp"

#endif  //ADF_PPC_HPP