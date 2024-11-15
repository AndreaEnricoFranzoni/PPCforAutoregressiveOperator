#include "strategy_cv.hpp"


//Augmenting window strategy division in train and validation
template<CV_STRAT cv_strat>
void
cv_strategy<cv_strat>::train_validation_set_strategy(int min_dim_ts, int max_dim_ts, CV_STRAT_T<CV_STRAT::AUGMENTING_WINDOW>)
{
  
  //defining the augmenting moving window for the train and validation sets 
  //how many, from the start, contiguos time instants: taking the ceil of half of the time instants as the smallest one
  //int min_dim_ts = static_cast<int>(std::ceil(static_cast<double>(n)/static_cast<double>(2)));
  //the biggest is all the time instants but one
  //int max_dim_ts = n;

  std::size_t min_size_ts = static_cast<std::size_t>(min_dim_ts);
  std::size_t max_size_ts = static_cast<std::size_t>(max_dim_ts);
  
  //train set is: from the beginning up to a time instant
  //validation set is: the next time instant
  m_strategy.reserve(max_size_ts - min_size_ts);
  for(std::size_t i = min_size_ts; i < max_size_ts; ++i)
  { 
    //train set: i: Eigen takes the first p left cols starting from 1
    std::vector<int> train_set;
    train_set.reserve(1);
    train_set.emplace_back(i);
    
    //validation set: i: Eigen takes the i-th column starting from 0
    std::vector<int> validation_set;
    validation_set.reserve(1);
    validation_set.emplace_back(i);
    
    m_strategy.emplace_back(std::make_pair(train_set,validation_set));
  }
}



//Augmenting window strategy test and validation set
//the first element of the pair is the test set, the second one is validation set
template<CV_STRAT cv_strat>
train_valid_set_t
cv_strategy<cv_strat>::train_validation_set(const KO_Traits::StoringMatrix &data, const iter_cv_t &strat, CV_STRAT_T<CV_STRAT::AUGMENTING_WINDOW>)
const
{
  return std::make_pair( data.leftCols(strat.first.front()), data.col(strat.second.front()) );
}