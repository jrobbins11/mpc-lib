/*
based heavily on osqp-eigen example: 
https://github.com/robotology/osqp-eigen/tree/master/example/src

and on Borrelli "Predictive Control for Linear and Hybrid Systems" 11.3.2
*/

/* 
TO DO: implement some sort of fault behavior if OsqpEigen returns false
for any of its methods
*/

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "OsqpEigen/OsqpEigen.h"
#include <stdexcept>
#include <vector>
#include <iostream>
#include <limits>

class MpcController
{
  private:
    // utilities
    void getTripletsForMatrix(const Eigen::MatrixXd &mat, 
      std::vector<Eigen::Triplet<double>> &tripvec,
      int rowOffset, int colOffset);

    // infinity
    const double inf = std::numeric_limits<double>::max();

  protected:

    // number of states, inputs, and inequality constraints
    int n_states;
    int n_inputs;
    int n_xineq;
    int n_uineq;

    // prediction horizon
    int n_horizon; 

    // initial condition
    Eigen::VectorXd x0;

    // reference
    Eigen::MatrixXd x_ref;

    // dynamics matrices
    Eigen::MatrixXd A_dyn;
    Eigen::MatrixXd B_dyn;

    // cost function matrices
    Eigen::MatrixXd Q_cost;
    Eigen::MatrixXd R_cost;
    Eigen::MatrixXd P_cost; // terminal cost

    // inequality matrices
    // TO DO: rewrite this in terms of upper and lower bounds so that matrices will be smaller
    Eigen::MatrixXd Ax_ineq;
    Eigen::VectorXd bx_ineq_low;
    Eigen::VectorXd bx_ineq_up;
    Eigen::MatrixXd Au_ineq;
    Eigen::VectorXd bu_ineq_low;
    Eigen::VectorXd bu_ineq_up;

    // optimization problem matrices
    Eigen::SparseMatrix<double> H; // hessian
    Eigen::VectorXd f; // gradient
    Eigen::SparseMatrix<double> A; // constraints
    Eigen::VectorXd b_low; // lower bound
    Eigen::VectorXd b_up; // upper bound

    // solver variables
    OsqpEigen::Solver solver;
    Eigen::VectorXd solution;

    // setup optimization problem
    void makeInequalityMatrices();
    void makeCostMatrices();

    // solver functions
    bool configureSolver();
    bool solveOptimizationProblem();

  public:

    // constructor
    MpcController(const Eigen::VectorXd &x0, const Eigen::MatrixXd &x_ref, 
      const Eigen::MatrixXd &A_dyn, const Eigen::MatrixXd &B_dyn,
      const Eigen::MatrixXd &Q_cost, const Eigen::MatrixXd &R_cost, 
      const Eigen::MatrixXd &P_cost, const Eigen::MatrixXd &Ax_ineq, 
      const Eigen::VectorXd &bx_ineq_low, const Eigen::VectorXd &bx_ineq_up,
      const Eigen::MatrixXd &Au_ineq, const Eigen::VectorXd &bu_ineq_low,
      const Eigen::VectorXd &bu_ineq_up, int n_horizon);

    // control method
    Eigen::VectorXd control(const Eigen::VectorXd &x, const Eigen::MatrixXd &x_ref);

    // print method
    void printOptimizationProblem();

};