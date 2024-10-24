#ifndef CV_CRTP_ALPHA_K_PPC_HPP
#define CV_CRTP_ALPHA_K_PPC_HPP

#include "CV.hpp"
#include "CV_alpha_k.hpp"


template< CV_STRAT cv_strat, CV_ERR_EVAL err_eval, K_IMP k_imp, VALID_ERR_RET valid_err_ret >  
class CV_alpha_k : public CV_base< CV_alpha_k<cv_strat,err_eval,k_imp,valid_err_ret>, cv_strat, err_eval, k_imp, valid_err_ret >
{
private:
  std::vector<double> m_alphas;                   //vector of alphas
  std::vector<int> m_k_s;                         //vector of k_s
  double m_alpha_best;                            //best alpha value
  int m_k_best;                                   //best k value
  double m_best_valid_error;                      //validation error for the best pair
  std::map<double,int> m_best_pairs;              //each alphas with its best k
  valid_err_cv_2_t m_valid_errors;                //to save the error for each pair alpha-k tried
  std::vector<double> m_valid_errors_best_pairs;  //the error in the position i-th is relative to the pair alpha-k in the map, where alpha is the i-th in m_alphas, and k is the one that gives the best validation error given the alpha
  double m_toll;
  
  pred_func_t<K_IMP::YES> m_pred_f;               //function to make predictions
  
public:
  //constructor
  template<typename STOR_OBJ,typename STRATEGY>
  CV_alpha_k(STOR_OBJ&& Data,
             STRATEGY && strategy,
             const std::vector<double> &alphas,
             const std::vector<int> &k_s,
             double toll,
             const pred_func_t<K_IMP::YES> & pred_f)
    : CV_base<CV_alpha_k,cv_strat,err_eval,k_imp,valid_err_ret>(std::move(Data),std::move(strategy)), 
      m_alphas(alphas),
      m_k_s(k_s),
      m_toll(toll),
      m_pred_f(pred_f)
      {
        
         /*
          * std::cout << "CV alpha-k" << std::endl;
          std::cout << "Data" << std::endl;
          std::cout << this->Data() << std::endl;
          std::cout << "tol: " << m_toll << std::endl;
          std::cout << "alphas:" << std::endl;
          for(std::size_t i = 0; i < m_alphas.size(); ++i){std::cout<<m_alphas[i]<<std::endl;}
          std::cout << "k_s: " << std::endl;
          for(std::size_t i = 0; i < m_k_s.size(); ++i){std::cout<<m_k_s[i]<<std::endl;}
          */
         
      }
  
  
  //Getters
  /*!
   * Getter for m_alpha_best
   */
  inline double alpha_best() const {return m_alpha_best;};
  
  /*!
   * Getter for m_k_best
   */
  inline double k_best() const {return m_k_best;};
  
  /*!
   * Getter for m_valid_errors_best_pairs
   */
  inline std::vector<double> valid_errors_best_pairs() const {return m_valid_errors_best_pairs;};
  
  /*!
   * Getter for m_valid_errors
   */
  inline std::vector<std::vector<double>> valid_errors() const {return m_valid_errors;};
  
  
  /*!
   * Getter for m_best_valid_error
   */
  inline double best_valid_error() const {return m_best_valid_error;};
  
  
  
  inline 
  void 
  best_param_search()
  {
    //preparing the containers for the errors
    if constexpr(valid_err_ret == VALID_ERR_RET::YES_err)
    {
      m_valid_errors.reserve(m_alphas.size());  //for each alpha: valid error for each k
    }
    
    m_valid_errors_best_pairs.reserve(m_alphas.size()); //for each alpha: valid error only for the best k
    

    for(const auto & alpha : m_alphas)
    {
      //alpha fixed: doing CV on k
      CV_k<cv_strat,err_eval,k_imp,valid_err_ret> cv(std::move(this->Data()),std::move(this->strategy()),m_k_s,m_toll,alpha,m_pred_f);
      cv.best_param_search();
      
      //best k given the alpha
      m_best_pairs.insert(std::make_pair(alpha,cv.param_best()));
      //saving the validation error for the best pair
      m_valid_errors_best_pairs.emplace_back(cv.best_valid_error());
      
      if constexpr(valid_err_ret == VALID_ERR_RET::YES_err)
      {
        //saving the validation error for each pair (further inspection)
        m_valid_errors.emplace_back(cv.valid_errors());
      }
      
    }
    
    //best validation error
    auto min_err = std::min_element(m_valid_errors_best_pairs.begin(),m_valid_errors_best_pairs.end());
    m_best_valid_error = *min_err;
    
    //best alpha
    m_alpha_best = m_alphas[std::distance(m_valid_errors_best_pairs.begin(),min_err)];
    
    //best k
    m_k_best = m_best_pairs.find(m_alpha_best)->second;
  }
  
};
  

#endif  //CV_CRTP_ALPHA_K_PPC_HPP