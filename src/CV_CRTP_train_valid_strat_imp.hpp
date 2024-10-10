#include "CV_CRTP.hpp"







//Augmenting window strategy division in train and validation
template<class D, DEF_PARAMS_PPC::cv_strat_type cv_strat, DEF_PARAMS_PPC::cv_err_eval_type err_eval>
cv_strategy_t
CV_PPC::CV_KO_PPC_CRTP<D,cv_strat,err_eval>::train_validation_set_strategy(const KO_Traits::StoringMatrix &dataset,CV_STRAT_T<DEF_PARAMS_PPC::cv_strat_type::AUGMENTING_WINDOW>)
const
{
  
  //defining the augmenting moving window for the train and validation sets 
  //how many, from the start, contiguos time instants: taking the ceil of half of the time instants as the smallest one
  int min_dim_ts = static_cast<int>(std::ceil(static_cast<double>(dataset.cols())/static_cast<double>(2)));
  //the biggest is all the time instants but one
  int max_dim_ts = static_cast<int>(dataset.cols());
  
  //train set is: from the beginning up to a time instant
  //validation set is: the next time instant
  cv_strategy_t strategy;
  strategy.reserve(max_dim_ts - min_dim_ts);
  for(std::size_t i = static_cast<std::size_t>(min_dim_ts); i < static_cast<std::size_t>(max_dim_ts); ++i)
  { 
    //train set: i: Eigen takes the first p left cols starting from 1
    std::vector<int> train_set;
    train_set.reserve(1);
    train_set.emplace_back(i);
    
    //validation set: i: Eigen takes the i-th column starting from 0
    std::vector<int> validation_set;
    validation_set.reserve(1);
    validation_set.emplace_back(i);
    
    strategy.emplace_back(std::make_pair(train_set,validation_set));
  }
  
  return strategy;
}








//Augmenting window strategy test and validation set
//the first element of the pair is the test set, the second one is validation set
template<class D, DEF_PARAMS_PPC::cv_strat_type cv_strat, DEF_PARAMS_PPC::cv_err_eval_type err_eval>
train_valid_set_t
CV_PPC::CV_KO_PPC_CRTP<D,cv_strat,err_eval>::train_validation_set(const std::pair<std::vector<int>,std::vector<int>> &strat,CV_STRAT_T<DEF_PARAMS_PPC::cv_strat_type::AUGMENTING_WINDOW>)
const
{
  return std::make_pair( this->Data().leftCols(strat.first.front()), this->Data().col(strat.second.front()) );
}
  