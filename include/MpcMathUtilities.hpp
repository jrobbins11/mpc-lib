#ifndef _MPCMATHUTILITIES_HPP_
#define _MPCMATHUTILITIES_HPP_

#include <cmath>
#include "Eigen/Dense"
#include "Eigen/Sparse"

// factorial
unsigned int factorial(unsigned int x)
{
    return (x==0) ? 1 : x*factorial(x-1);
}

// matrix power
Eigen::MatrixXd pow(const Eigen::MatrixXd &M, unsigned int N) 
{
  // check dimensions
  if (M.rows() != M.cols())
    throw std::invalid_argument("M is not square");

  // compute matrix power
  Eigen::MatrixXd Mpow = Eigen::MatrixXd::Identity(M.rows(), M.cols()); // init
  for (unsigned int i=0; i<N; i++)
  {
    Mpow *= M;
  }
  return Mpow;
}

#endif