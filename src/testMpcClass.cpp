#include <iostream>
#include <limits>
#include <ctime>
#include "MpcController.hpp"

Eigen::VectorXd stepDynamics(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B,
  const Eigen::VectorXd &x, const Eigen::VectorXd &u)
{
    return A*x + B*u;
}

int main()
{
    // constants
    const double inf = std::numeric_limits<double>::max();

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

    // state constraint matrices: x1 in [-1, 1], x2 in [-1, 1]
    Eigen::Matrix<double, 2, 2> I2 = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 2, 2> Ax_ineq;
    Eigen::Vector<double, 2> bx_ineq_low;
    Eigen::Vector<double, 2> bx_ineq_up;
    Ax_ineq << I2;
    bx_ineq_low << -1, -1;
    bx_ineq_up << 1, 1;

    // terminal state constraint matrices: same as state constraint matrices here
    Eigen::Matrix<double, 2, 2> Ax_term_ineq;
    Eigen::Vector<double, 2> bx_term_ineq_low;
    Eigen::Vector<double, 2> bx_term_ineq_up;
    Ax_term_ineq << I2;
    bx_term_ineq_low << -1, -1;
    bx_term_ineq_up << 1, 1;

    // state constraint softening
    Eigen::DiagonalMatrix<double, 2> Qx_constraint_cost;
    Qx_constraint_cost.diagonal() << 1e6, 1e6;

    Eigen::DiagonalMatrix<double, 2> Qxterm_constraint_cost; // terminal constraint
    Qxterm_constraint_cost.diagonal() << 1e6, 1e6;

    Eigen::DiagonalMatrix<double, 1> Qu_constraint_cost;
    Qu_constraint_cost.diagonal() << 1e6;

    // input constraint matrices: u in [-1, 1]
    Eigen::Matrix<double, 1, 1> Au_ineq;
    Eigen::Vector<double, 1> bu_ineq_low;
    Eigen::Vector<double, 1> bu_ineq_up;
    Au_ineq << 1;
    bu_ineq_low << -1;
    bu_ineq_up << 1;

    // prediction horizon
    int n_horizon = 10;

    // loop time of controller
    const double t_loop = 0.1;

    // initial condition
    Eigen::Vector<double, 2> x0;
    x0 << 1, 0;

    // references, allowed to vary over horizon
    Eigen::MatrixXd x_ref = Eigen::MatrixXd::Zero(2, n_horizon);

    // LTV setup
    std::vector<Eigen::MatrixXd> A_dyn_vec;
    std::vector<Eigen::MatrixXd> B_dyn_vec;
    for (int i=0; i<n_horizon; i++)
    {
      A_dyn_vec.push_back(A_dyn);
      B_dyn_vec.push_back(B_dyn);
    }

    /*
    // create MPC object - no constraint softening
    MpcController MPC(A_dyn, B_dyn, 
      Q_cost, R_cost, P_cost, 
      Ax_ineq, bx_ineq_low, bx_ineq_up, 
      Ax_term_ineq, bx_term_ineq_low, bx_term_ineq_up, 
      Au_ineq, bu_ineq_low, bu_ineq_up, 
      n_horizon, t_loop);
    */

    /*
    // create MPC object - with softened state constraints
    MpcController MPC(A_dyn, B_dyn, 
      Q_cost, R_cost, P_cost, 
      Qx_constraint_cost, Qxterm_constraint_cost,
      Ax_ineq, bx_ineq_low, bx_ineq_up, 
      Ax_term_ineq, bx_term_ineq_low, bx_term_ineq_up, 
      Au_ineq, bu_ineq_low, bu_ineq_up, 
      n_horizon, t_loop);
    */

    /*
    // create MPC object - with softened state and input constraints
    MpcController MPC(A_dyn, B_dyn, 
      Q_cost, R_cost, P_cost, 
      Qx_constraint_cost, Qxterm_constraint_cost, Qu_constraint_cost,
      Ax_ineq, bx_ineq_low, bx_ineq_up, 
      Ax_term_ineq, bx_term_ineq_low, bx_term_ineq_up, 
      Au_ineq, bu_ineq_low, bu_ineq_up, 
      n_horizon, t_loop);
    */

    // build MPC object from default constructor
    MpcController MPC;
    //MPC.setDynMatrices(A_dyn, B_dyn);
    MPC.setDynMatrices(A_dyn_vec, B_dyn_vec); // LTV
    MPC.setStageCost(Q_cost, R_cost);
    MPC.setMpcHorizon(n_horizon);
    MPC.setMaxExecutionTime(t_loop);
    MPC.buildController();

    // display optimization problem matrices
    //MPC.printOptimizationProblem();

    // declare control input variable
    Eigen::Vector<double, 1> u; 

    // time execution
    clock_t timer;
    float avg_time, time;
    float tot_time = 0;
    float max_time = 0;

    // run MPC controller in loop and print control inputs to screen
    int n_cycles = 19; // numbe of cycles 
    Eigen::Vector<double, 2> x = x0;
    for (int i=0; i<n_cycles; i++)
    {
     
      // start timer
      timer = clock();

      // compute control
      // u = MPC.control(x, x_ref);
      u = MPC.control(x, x_ref, A_dyn_vec, B_dyn_vec); // LTV

      // log execution time
      timer = clock() - timer;
      time = (float) timer/CLOCKS_PER_SEC;
      tot_time += time;
      if (time > max_time)
        max_time = time;

      // print control output to console
      std::cout << "u = " << u << std::endl;

      // step dynamics
      x = stepDynamics(A_dyn, B_dyn, x, u);
    }
    
    // display execution time
    std::cout << "Average MPC execution time = " << 
        tot_time/n_cycles << " sec" << std::endl;
    std::cout << "Max MPC execution time = " <<
        max_time << " sec" << std::endl;

    return 0;
}