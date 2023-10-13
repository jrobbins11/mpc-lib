#include <iostream>
#include "MpcController.hpp"

Eigen::VectorXd stepDynamics(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B,
  const Eigen::VectorXd &x, const Eigen::VectorXd &u)
{
    return A*x + B*u;
}

int main()
{
    // dynamics matrices
    Eigen::Matrix<double, 2, 2> A_dyn;
    Eigen::Matrix<double, 2, 1> B_dyn;
    A_dyn << 1, 1, 0, 1;
    B_dyn << 0, 1;

    // cost matrices
    Eigen::DiagonalMatrix<double, 2> Q_cost;
    Eigen::DiagonalMatrix<double, 1> R_cost;
    Eigen::DiagonalMatrix<double, 2> P_cost;
    Q_cost.diagonal() << 1, 1;
    R_cost.diagonal() << 0.1;
    P_cost.diagonal() << 1, 1;

    // state constraint matrices: x in [-1, 1]
    Eigen::Matrix<double, 2, 2> I2 = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 4, 2> Ax_ineq;
    Eigen::Vector<double, 4> bx_ineq;
    Ax_ineq << I2,
              -I2; // writing this way for readability
    bx_ineq << 1, 1, 1, 1;

    // input constraint matrices: u in [-1, 1]
    Eigen::Matrix<double, 2, 1> Au_ineq;
    Eigen::Vector<double, 2> bu_ineq;
    Au_ineq << 1, -1;
    bu_ineq << 1, 1;

    // prediction horizon
    int n_horizon = 10;

    // initial condition
    Eigen::Vector<double, 2> x0;
    x0 << 1, 0;

    // references, allowed to vary over horizon
    Eigen::MatrixXd x_ref = Eigen::MatrixXd::Zero(2, n_horizon);

    // create MPC object
    MpcController MPC(x0, x_ref, A_dyn, B_dyn, Q_cost, R_cost, P_cost, 
      Ax_ineq, bx_ineq, Au_ineq, bu_ineq, n_horizon);

    // display optimization problem matrices
    //MPC.printOptimizationProblem();

    // declare control input variable
    Eigen::Vector<double, 1> u; 

    // run MPC controller in loop and print control inputs to screen
    Eigen::Vector<double, 2> x = x0;
    for (int i=0; i<19; i++)
    {
      // compute control and print to console
      u = MPC.control(x, x_ref);
      std::cout << "u = " << u << std::endl;

      // step dynamics
      x = stepDynamics(A_dyn, B_dyn, x, u);
    }
    
    return 0;
}