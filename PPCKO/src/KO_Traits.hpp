#ifndef KO_TRAITS_HPP
#define KO_TRAITS_HPP

#include <Eigen/Dense>


// Defining common types for homogeneity
struct KO_Traits
{
public:
  
  // Data are stored in dynamic matrices (sizes are known only at compile time, and easily very larges)
  // of doubles: the only conversion from R types is Numeric -> double. some precision is lost. But double are
  // used since time series are of real numbers.
  using StoringMatrix = Eigen::MatrixXd;
  
  using StoringVector = Eigen::VectorXd;    //col-wise
  
  using StoringArray = Eigen::ArrayXd;
  
};

#endif /*KO_TRAITS_HPP*/