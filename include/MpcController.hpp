/*
based heavily on osqp-eigen example: 
https://github.com/robotology/osqp-eigen/tree/master/example/src

and on Borrelli "Predictive Control for Linear and Hybrid Systems" 11.3.2
*/

/* 
TO DO: implement some sort of fault behavior if OsqpEigen returns false
for any of its methods
*/

#ifndef _MPCCONTROLLER_HPP_
#define _MPCCONTROLLER_HPP_

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
    int n_states = 0;
    int n_inputs = 0;
    int n_xineq = 0;
    int n_uineq = 0;
    int n_xtermineq = 0;

    // prediction horizon
    int n_horizon = 0; 

    // loop time
    double T_sec = inf;

    // flags
    bool terminalCostSpecified = false;
    bool softenedStateConstraints = false;
    bool softenedInputConstraints = false;
    bool stateConstraintSpecified = false;
    bool inputConstraintSpecified = false;
    bool terminalStateConstraintSpecified = false;
    bool terminalStateConstraintSofteningCostSpecified = false;
    bool LTV = false;

    // initial condition
    Eigen::VectorXd x0 = Eigen::VectorXd::Zero(0);

    // reference
    Eigen::MatrixXd x_ref = Eigen::MatrixXd::Zero(0,0);

    // dynamics matrices
    Eigen::MatrixXd A_dyn = Eigen::MatrixXd::Zero(0,0);
    Eigen::MatrixXd B_dyn = Eigen::MatrixXd::Zero(0,0);
    std::vector<Eigen::MatrixXd> A_dyn_vec;
    std::vector<Eigen::MatrixXd> B_dyn_vec;

    // cost function matrices
    Eigen::MatrixXd Q_cost = Eigen::MatrixXd::Zero(0,0);
    Eigen::MatrixXd R_cost = Eigen::MatrixXd::Zero(0,0);
    Eigen::MatrixXd P_cost = Eigen::MatrixXd::Zero(0,0); // terminal cost

    // softened constraint costs
    // J_k = (1/2)*(s_k^T) Q (s_k)
    // where bx_low <= Ax*x_k - sx_k <= bx_high
    // where bu_low <= Au*u_k - su_k <= bu_high
    Eigen::MatrixXd Qx_constraint_cost = Eigen::MatrixXd::Zero(0,0);
    Eigen::MatrixXd Qxterm_constraint_cost = Eigen::MatrixXd::Zero(0,0);
    Eigen::MatrixXd Qu_constraint_cost = Eigen::MatrixXd::Zero(0,0);

    // state constraints
    Eigen::MatrixXd Ax_ineq = Eigen::MatrixXd::Zero(0,0);
    Eigen::VectorXd bx_ineq_low = Eigen::VectorXd::Zero(0);
    Eigen::VectorXd bx_ineq_up = Eigen::VectorXd::Zero(0);

    // terminal state constraints
    Eigen::MatrixXd Ax_term_ineq = Eigen::MatrixXd::Zero(0,0);
    Eigen::VectorXd bx_term_ineq_low = Eigen::VectorXd::Zero(0);
    Eigen::VectorXd bx_term_ineq_up = Eigen::VectorXd::Zero(0);

    // input constraints
    Eigen::MatrixXd Au_ineq = Eigen::MatrixXd::Zero(0,0);
    Eigen::VectorXd bu_ineq_low = Eigen::VectorXd::Zero(0);
    Eigen::VectorXd bu_ineq_up = Eigen::VectorXd::Zero(0);

    // optimization problem matrices
    Eigen::SparseMatrix<double> H; // hessian
    Eigen::VectorXd f; // gradient
    Eigen::SparseMatrix<double> A; // constraints
    Eigen::VectorXd b_low; // lower bound
    Eigen::VectorXd b_up; // upper bound

    // solver variables
    OsqpEigen::Solver solver;
    Eigen::VectorXd solution = Eigen::VectorXd::Zero(0);

    // setup optimization problem
    void makeInequalityMatrices();
    void makeCostMatrices();

    // solver functions
    bool configureSolver();
    bool solveOptimizationProblem();

    // utilities
    void validityCheckDimensions();

  public:

    // constructor - default
    MpcController() = default;

    // constructor - no constraint softening
    MpcController(
      const Eigen::MatrixXd &A_dyn, 
      const Eigen::MatrixXd &B_dyn,
      const Eigen::MatrixXd &Q_cost, 
      const Eigen::MatrixXd &R_cost, 
      const Eigen::MatrixXd &P_cost, 
      const Eigen::MatrixXd &Ax_ineq, 
      const Eigen::VectorXd &bx_ineq_low, 
      const Eigen::VectorXd &bx_ineq_up,
      const Eigen::MatrixXd &Ax_term_ineq, 
      const Eigen::VectorXd &bx_term_ineq_low, 
      const Eigen::VectorXd &bx_term_ineq_up,
      const Eigen::MatrixXd &Au_ineq, 
      const Eigen::VectorXd &bu_ineq_low,
      const Eigen::VectorXd &bu_ineq_up, 
      int n_horizon, 
      double T_loop_sec);

    // constructor - softened state constraints
    MpcController(
      const Eigen::MatrixXd &A_dyn, 
      const Eigen::MatrixXd &B_dyn,
      const Eigen::MatrixXd &Q_cost, 
      const Eigen::MatrixXd &R_cost, 
      const Eigen::MatrixXd &P_cost, 
      const Eigen::MatrixXd &Qx_constraint_cost,
      const Eigen::MatrixXd &Qxterm_constraint_cost,
      const Eigen::MatrixXd &Ax_ineq, 
      const Eigen::VectorXd &bx_ineq_low, 
      const Eigen::VectorXd &bx_ineq_up, 
      const Eigen::MatrixXd &Ax_term_ineq, 
      const Eigen::VectorXd &bx_term_ineq_low, 
      const Eigen::VectorXd &bx_term_ineq_up,
      const Eigen::MatrixXd &Au_ineq, 
      const Eigen::VectorXd &bu_ineq_low, 
      const Eigen::VectorXd &bu_ineq_up,
      int n_horizon, 
      double T_loop_sec);

    // constructor - softened state and input constraints
    MpcController(
      const Eigen::MatrixXd &A_dyn, 
      const Eigen::MatrixXd &B_dyn,
      const Eigen::MatrixXd &Q_cost, 
      const Eigen::MatrixXd &R_cost, 
      const Eigen::MatrixXd &P_cost, 
      const Eigen::MatrixXd &Qx_constraint_cost, 
      const Eigen::MatrixXd &Qxterm_constraint_cost,
      const Eigen::MatrixXd &Qu_constraint_cost,
      const Eigen::MatrixXd &Ax_ineq, 
      const Eigen::VectorXd &bx_ineq_low, 
      const Eigen::VectorXd &bx_ineq_up, 
      const Eigen::MatrixXd &Ax_term_ineq, 
      const Eigen::VectorXd &bx_term_ineq_low, 
      const Eigen::VectorXd &bx_term_ineq_up,
      const Eigen::MatrixXd &Au_ineq, 
      const Eigen::VectorXd &bu_ineq_low, 
      const Eigen::VectorXd &bu_ineq_up,
      int n_horizon, 
      double T_loop_sec);

    // problem setup methods
    void setDynMatrices(const Eigen::MatrixXd &A_dyn, const Eigen::MatrixXd &B_dyn);
    void setDynMatrices(const std::vector<Eigen::MatrixXd> &A_dyn_vec, const std::vector<Eigen::MatrixXd> &B_dyn_vec);
    void setStageCost(const Eigen::MatrixXd &Q_cost, const Eigen::MatrixXd &R_cost);
    void setTerminalCost(const Eigen::MatrixXd &P_cost);
    void setStateConstraints(const Eigen::MatrixXd &Ax_ineq, 
      const Eigen::MatrixXd &bx_ineq_low, const Eigen::MatrixXd &bx_ineq_up);
    void setInputConstraints(const Eigen::MatrixXd &Au_ineq, 
      const Eigen::MatrixXd &bu_ineq_low, const Eigen::MatrixXd &bu_ineq_up);
    void setTerminalStateConstraints(const Eigen::MatrixXd &Ax_term_ineq, 
      const Eigen::MatrixXd &bx_term_ineq_low, const Eigen::MatrixXd &bx_term_ineq_up);
    void setStateConstraintSlackVarsCost(const Eigen::MatrixXd &Qx_constraint_cost);
    void setInputConstraintSlackVarsCost(const Eigen::MatrixXd &Qu_constraint_cost);
    void setStateTerminalConstraintSlackVarsCost(const Eigen::MatrixXd &Qxterm_constraint_cost);
    void setMpcHorizon(int n_horizon);
    void setMaxExecutionTime(double T_loop_sec);
    void buildController();

    // control method
    Eigen::VectorXd control(const Eigen::VectorXd &x, const Eigen::MatrixXd &x_ref);
    Eigen::VectorXd control(const Eigen::VectorXd &x, const Eigen::MatrixXd &x_ref,
      const std::vector<Eigen::MatrixXd> &A_dyn_vec, const std::vector<Eigen::MatrixXd> &B_dyn_vec);

    // print methods
    void printOptimizationProblem();
    void printSolution();

};

#endif