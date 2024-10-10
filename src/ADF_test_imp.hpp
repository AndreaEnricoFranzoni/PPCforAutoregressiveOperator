#include "ADF_test.hpp"



//function to evaluate the one time step differences of the time series
template<class LAG_policy>
std::vector<double>
PPC_util::adf<LAG_policy>::one_step_diff(const KO_Traits::StoringVector &ts)
const
{
  std::vector<double> diff;
  diff.resize(m_tot_time_instants-1);
  std::adjacent_difference(ts.cbegin()+1,ts.cend(),diff.begin());
  diff[0] = ts[1] - ts[0];
  
  return diff;
}



//function to embed the time series x into a low-dimensional Euclidean space
template<class LAG_policy>
KO_Traits::StoringMatrix
PPC_util::adf<LAG_policy>::embed(const KO_Traits::StoringVector &ts, int dimension)
const
{ 
  //evaluate one step differences for the time serie
  std::vector<double> y(this->one_step_diff(ts));

  int n = y.size();
  int m = n - dimension + 1;    //number of rows of embedded ts
  
  KO_Traits::StoringMatrix embedded(m,dimension);
  
  //embedding
  std::vector<std::size_t> start_ind;
  start_ind.resize(static_cast<std::size_t>(m));
  std::iota(start_ind.begin(),start_ind.end(),static_cast<std::size_t>(0));
  
  for(std::size_t i = 0; i < static_cast<std::size_t>(dimension); ++i)
  {
    std::size_t gain = static_cast<std::size_t>(dimension) - i - 1;
    
    std::transform(start_ind.begin(),start_ind.end(),
                   embedded.col(i).begin(),
                   [&y,&gain](std::size_t el){ return y[el+gain];});
    
  }
  
  start_ind.clear();
  
  return embedded;
}



//evaluation of the test statistics depending on the lag
template<class LAG_policy>
double
PPC_util::adf<LAG_policy>::statistic_eval(const KO_Traits::StoringVector &ts)
const
{
  //embedding of the time serie in an Euclidean space
  KO_Traits::StoringMatrix z = this->embed(ts,m_k_used);
  
  //the design matrix depends on the time lag: Lag_policy will take care of it
  ADF_util::adf_stat<LAG_policy> statistic_adf;
  return statistic_adf(ts,z);
}
  


//evaluate the pvalue for the test on a single time serie
template<class LAG_policy>
double
PPC_util::adf<LAG_policy>::p_value_eval(const KO_Traits::StoringVector &ts, const std::vector<double> &tableipl, int i) 
const
{ 
  //evaluation of the test statistic
  double stat = this->statistic_eval(ts);
  
  //if statistic too extreme the pvalue is put to 0 or 1
  if(stat<tableipl.front()){return 0.0;}
  if(stat>tableipl.back()){return 1.0;}
  
  //pvalue evaluation
  PPC_util::integrand_interp p_val_int{tableipl};
  double pval = p_val_int(stat,ADF_util::tablep);
  
  return pval;
}



//doing the test for each temporal series
template<class LAG_policy>
void
PPC_util::adf<LAG_policy>::test()
{
  
   //preparing for the p_value evaluation
   m_p_values.reserve(m_x.rows());
   std::vector<double> tableipl = ADF_util::tableipl(static_cast<double>(m_tot_time_instants-1));
   
   //computing the pvalues for each time serie
   for (int i = 0; i < m_x.rows(); ++i) 
   {
      m_p_values.emplace_back( this->p_value_eval(m_x.row(i),tableipl,i+1)  );
   }
   
}