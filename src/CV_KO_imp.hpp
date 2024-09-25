#include "CV_KO.hpp"


template<typename T>
double
CV_PPC::CV_KO<T>::single_cv(T param,
                            int dim_train_set) 
const
{   
  
  if constexpr(std::is_same<T,double>::value)
  { 
    //a single ko with a given alpha parameter on a specific training set
    //prediction with mean already added on non-normalized real values
    return this->ef()( this->singleCV()(this->X().leftCols(dim_train_set), m_threshold_ppc, param, m_k).array()     - m_X.col(dim_train_set).array() );
  }
  else if constexpr(std::is_same<T,int>::value)
  { 
    //a single ko with a given alpha parameter on a specific training set
    //prediction with mean already added on non-normalized real values
    return this->ef()( this->singleCV()(this->X().leftCols(dim_train_set), m_threshold_ppc, m_alpha, param).array() - m_X.col(dim_train_set).array() );
  }
}


//mean of all the mse, alpha fixed, given all the exploited training set
template<typename T>
double
CV_PPC::CV_KO<T>::moving_window_cv(T param,
                                   const std::vector<int> & t_i) 
const
{   

  //given a param, evaluating all the ef, and summing them
  double errors = std::transform_reduce(t_i.cbegin(),
                                        t_i.cend(),
                                        0.0,
                                        std::plus{},
                                        [this,&param](int const &i){return this->single_cv(param,i);});
    
    //returning the mean of the ef
    return errors/static_cast<double>(t_i.size()); 
}



template<typename T>
void
CV_PPC::CV_KO<T>::best_param()
{   
  //defining the moving window for the train and validation sets 
  //how many, from the start, contiguos time instants: taking the ceil of half of the time instants as the smallest one
  int min_dim_ts = static_cast<int>(std::ceil(static_cast<double>(this->X().cols())/static_cast<double>(2)));
  //the biggest is all the time instants but one
  int max_dim_ts = static_cast<int>(this->X().cols());
  
  //putting the time instants that indicate the end of the training set in this vector
  std::vector<int> time_instants_cv;
  time_instants_cv.resize(max_dim_ts-min_dim_ts);
  std::iota(time_instants_cv.begin(),time_instants_cv.end(),min_dim_ts);
  
  
  
  if constexpr(std::is_same<T,double>::value)     //alpha
  { 
    //looping over all the parameters
    std::transform(m_params.cbegin(),
                   m_params.cend(),
                   m_errors.begin(),
                   [this,&time_instants_cv](T const &param_i){return this->moving_window_cv(param_i,time_instants_cv);});
    
  }
  else if constexpr(std::is_same<T,int>::value)   //k
  { 
    std::size_t counter_k = 0;
    double previous_error = static_cast<double>(0);
    
    //if adding another PPC does not improve too much the validation error: break
    for(const auto & el : m_params)
    {
      m_errors[counter_k] = this->moving_window_cv(el,time_instants_cv);
      if(std::abs(m_errors[counter_k] - previous_error) < DEF_PARAMS_PPC::toll_cv_k)     //already found a point in which you do not improve anymore
      {
        m_errors.resize(el);
        break;
      }
      else
      {
        previous_error = m_errors[counter_k];
        counter_k++;
      }
    }
  }
  
  //optimal param
  m_param_best = m_params[std::distance(m_errors.begin(),std::min_element(m_errors.begin(),m_errors.end()))];

  //free memory
  time_instants_cv.clear();
  m_params.clear();
  //m_errors.clear();
}