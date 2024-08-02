#include "CV_KO.hpp"



//mse between a given training and validation set, alpha fixed
double
  CV_PPC::CV_KO::single_cv(double alpha,int dim_train_set, const std::function<double(KO_Traits::StoringVector)> &lf) 
  const
  {   
    //a single ko with a given alpha parameter on a specific training set
    PPC::KO_NO_CV ko_single_cv(std::move(this->X().leftCols(dim_train_set)),alpha);
    ko_single_cv.solve();
    
    //prediction with mean already added on non-normalized real values
    return lf(ko_single_cv.prediction() - m_X.col(dim_train_set).array());
  }


//mean of all the mse, alpha fixed, given all the exploited training set
double
  CV_PPC::CV_KO::moving_window_cv(double alpha, const std::vector<int> & t_i, const std::function<double(KO_Traits::StoringVector)> &lf) 
  const
  {   
    //given alpha, evaluating all the mses
    double errors = std::transform_reduce(t_i.cbegin(),
                                          t_i.cend(),
                                          0.0,
                                          std::plus{},
                                          [this,&alpha,&lf](int const &i){return this->single_cv(alpha,i,lf);}
    );
    
    //returning the mean of the mses
    return errors/static_cast<double>(t_i.size()); 
  }


void
CV_PPC::CV_KO::best_param()
{   
  //defining the error function to be used to evaluate the error during the CV process (can be changed)
  std::function<double(KO_Traits::StoringVector)> ef = EF_PPC::mse<double>;
  
  //defining the moving window for the train and validation sets 
  //how many, from the start, contiguos time instants: taking the ceil of half of the time instants as the smallest one
  int min_dim_ts = static_cast<int>(std::ceil(static_cast<double>(this->X().cols())/static_cast<double>(2)));
  //the biggest is all the time instants but one
  int max_dim_ts = static_cast<int>(this->X().cols());
  
  //putting the time instants that indicate the end of the training set in this vector
  std::vector<int> time_instants_cv;
  time_instants_cv.resize(max_dim_ts-min_dim_ts);
  std::iota(time_instants_cv.begin(),time_instants_cv.end(),min_dim_ts);
  
  //looping over all the alphas
  std::transform(m_alphas.cbegin(),
                 m_alphas.cend(),
                 m_errors.begin(),
                 [this,&time_instants_cv,&ef](double const &alpha_i){return this->moving_window_cv(alpha_i,time_instants_cv,ef);}
  );
  
  /*
   for(std::size_t i = 0; i < m_alphas.size(); ++i)
   {
   std::cout << i << ": for alpha equal to " << m_alphas[i] << " the validation error is " << m_errors[i] << std::endl;
   }
   */
  
  //returning the optimal alpha
  m_alpha_best = m_alphas[std::distance(m_errors.begin(),std::min_element(m_errors.begin(),m_errors.end()))];
  std::cout << "Alpha opt is: " << m_alpha_best << std::endl;
  
  //free memory (these vector are useless from now on)
  time_instants_cv.clear();
  m_alphas.clear();
  //m_errors.clear();
}