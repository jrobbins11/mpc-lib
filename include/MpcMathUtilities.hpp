#ifndef _MPCMATHUTILITIES_HPP_
#define _MPCMATHUTILITIES_HPP_

#include <cmath>
#include "Eigen/Dense"

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

// continuous to discrete dynamic matrix conversion
void cont2DiscDynMatrices(const Eigen::MatrixXd &Ac, const Eigen::MatrixXd &Bc, double dT, 
  Eigen::MatrixXd &Ad, Eigen::MatrixXd &Bd)
{
  // Taylor series order
  const int N = 10; 

  // initialize discrete matrices
  Ad.resize(Ac.rows(), Ac.cols());
  Bd.resize(Bc.rows(), Bc.cols());
  Ad = Eigen::MatrixXd::Zero(Ac.rows(), Ac.cols());
  Bd = Eigen::MatrixXd::Zero(Bd.rows(), Bd.cols());

  // Taylor series approximation
  for (int k=0; k<=N; k++)
  {
    Ad = Ad + (1/((double) factorial(k)))*pow(Ac*dT,k);
    if (k >= 1)
      Bd = Bd + ((1/((double) factorial(k)))*pow(Ac,k-1)*pow(dT,k))*Bc;
  }
}

#endif