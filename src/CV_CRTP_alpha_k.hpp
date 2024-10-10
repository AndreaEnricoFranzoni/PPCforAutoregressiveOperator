#ifndef CV_CRTP_ALPHA_K_PPC_HPP
#define CV_CRTP_ALPHA_K_PPC_HPP

#include "CV_CRTP.hpp"
#include "CV_CRTP_alpha_k.hpp"

namespace CV_PPC
{
  
template<DEF_PARAMS_PPC::cv_strat_type cv_strat, DEF_PARAMS_PPC::cv_err_eval_type err_eval>  
class CV_KO_PPC_alpha_k_CRTP : public CV_KO_PPC_CRTP<CV_KO_PPC_alpha_k_CRTP<cv_strat,err_eval>, cv_strat, err_eval>
{
private:
  std::vector<double> m_alphas; //vector of alphas
  std::vector<int> m_k_s;       //vector of k_s
  double m_alpha_best;      //best alpha value
  int m_k_best;             //best k value
  double m_best_valid_error;      //validation error for the best pair
  std::map<double,int> m_best_pairs;                //each alphas with its best k
  std::vector<std::vector<double>> m_valid_errors;        //to save the error for each pair alpha-k
  std::vector<double> m_valid_errors_best_pairs;            //the error in the position i-th is relative to the pair alpha-k in the map, where alpha is the i-th in m_alphas, and k is the one that gives the best validation error given the alpha
  
  
  double m_toll;
  //parameters needed for KO
  double m_threshold_ppc;
  pred_func_t m_pred_f;             //function to make predictions
  
public:
  //constructor
  CV_KO_PPC_alpha_k_CRTP(KO_Traits::StoringMatrix&& Data,
                         const std::vector<double> &alphas,
                         const std::vector<int> &k_s,
                         double toll,
                         double threshold_ppc,
                         const pred_func_t & pred_f)
    : CV_KO_PPC_CRTP<CV_KO_PPC_alpha_k_CRTP, cv_strat, err_eval>(std::move(Data)), 
      m_alphas(alphas),
      m_k_s(k_s),
      m_toll(toll),
      m_threshold_ppc(threshold_ppc),
      m_pred_f(pred_f)
      {}
  
  
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
    m_valid_errors.reserve(m_alphas.size());  //for each alpha: valid error for each k
    m_valid_errors_best_pairs.reserve(m_alphas.size()); //for each alpha: valid error only for the best k
    
    int counter = 0;
    for(const auto & alpha : m_alphas)
    {
      //alpha fixed: doing CV on k
      CV_PPC::CV_KO_PPC_k_CRTP<cv_strat,err_eval> cv(std::move(this->Data()), m_k_s, m_toll, m_threshold_ppc, alpha, m_pred_f);
      cv.best_param_search();
      
      //best k given the alpha
      m_best_pairs.insert(std::make_pair(alpha,cv.param_best()));
      //saving the validation error for the best pair
      m_valid_errors_best_pairs.emplace_back(cv.best_valid_error());
      //saving the validation error for each pair (further inspection)
      m_valid_errors.emplace_back(cv.valid_errors());
      
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
  
}   //end namespace CV_PPC

#endif  //CV_CRTP_ALPHA_K_PPC_HPP