#ifndef SCORES_PPC_HPP
#define SCORES_PPC_HPP

#include <Eigen/Dense>
#include <vector>
#include <map>
#include <utility>

#include "mesh.hpp"
#include "L2_scalar_product.hpp"

#include <iostream>

namespace PPC
{
  
class scores
{
  
private:
Geometry::Domain1D m_domain;                                              //domain of the functional element
Geometry::Mesh1D m_grid_evaluations;                                      //values in the domain for which I already have the funciton value
Geometry::Mesh1D m_grid_integration;                                      //grid for doing the scalar product
std::map<std::size_t,std::vector<double>> m_integrand_evaluations;        //storing, for each direction, the evaluations
std::size_t m_k;                                                          //number of scores to be evaluated
std::vector<double> m_scores_evaluations;                                 //vector containing the value of the scores: length: k
  
  
public:
  scores(const Eigen::Matrix<double,Eigen::Dynamic,1> &f_n,
         const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &directions,
         double a,
         double b,
         int N
        )
    : 
    m_domain(a,b), 
    m_grid_evaluations(m_domain,static_cast<int>(f_n.size())-static_cast<int>(1)),
    m_grid_integration(m_domain,N),
    m_k(directions.cols())
  {
    
    //<f_n,a_i> for every i = {1,...,k}
    //preparing the discrete evaluations of the function f_n*a_i in order to integrate it
    for(std::size_t i = 0; i < m_k; ++i)
    {
      Eigen::Matrix<double,Eigen::Dynamic,1> prod = f_n.array()*directions.col(i).array();
      std::vector<double> prod_v(prod.data(), prod.data() + prod.size());       //using vectors to make it easier, since there are not other algebric matricial operations

      m_integrand_evaluations.insert(std::make_pair(i+static_cast<std::size_t>(1),prod_v));
    }
  
  }
  
  /*!
   * Getter for m_scores_evaluations
   */
  inline std::vector<double> scores_evaluations() const {return m_scores_evaluations;};
  
  //function to actually do the scalar product for each one of the directions
  void evaluating_scores();
  
};
  
} //end namespace PPC

#endif  //SCORES_PPC_HPP