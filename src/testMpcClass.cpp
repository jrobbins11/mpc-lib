#include <iostream>
#include <limits>
#include <ctime>
#include <cmath>
#include "MpcController.hpp"

Eigen::VectorXd stepDynamics(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B,
  const Eigen::VectorXd &x, const Eigen::VectorXd &u)
{
    return A*x + B*u;
}

// math utilities
unsigned int factorial(unsigned int x)
{
    return (x==0) ? 1 : x*factorial(x-1);
}

Eigen::MatrixXd pow(const Eigen::MatrixXd &M, unsigned int N) // matrix power
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

// Ackermann nonlinear dynamics
Eigen::Vector<double,6> ackermannDyn(const Eigen::Vector<double,6> &x, const Eigen::Vector<double,2> &u,
  double wheelBase, double dT)
{
  // state: x = [xp; yp; th; psi; v; psidot]
  // input: u = [a; psiddot]

  // unpack
  double xp, yp, th, psi, v, psidot, a, psiddot;

  xp = x[0];
  yp = x[1];
  th = x[2];
  psi = x[3];
  v = x[4];
  psidot = x[5];

  a = u[0];
  psiddot = u[1];

  // Ackermann kinematics
  double xpdot, ypdot, thdot;
  xpdot = v*cos(th);
  ypdot = v*sin(th);
  thdot = (1/wheelBase)*v*tan(psi);

  // assemble derivative
  Eigen::Vector<double, 6> xdot;
  xdot << xpdot, ypdot, thdot, psidot, a, psiddot;

  // discrete time dynamics
  return x + xdot*dT;

}

// linearized dynamics
void ackermannLinearizedDynMatrices(const Eigen::Vector<double,6> &x, double wheelBase, 
  Eigen::Matrix<double,6,6> &Ac, Eigen::Matrix<double,6,3> &Bc)
{
  // unpack
  double xp, yp, th, psi, v, psidot;

  xp = x[0];
  yp = x[1];
  th = x[2];
  psi = x[3];
  v = x[4];
  psidot = x[5];

  // A matrix
  Eigen::Matrix<double, 6, 6> A;
  A << 0, 0, -v*sin(th), 0, cos(th), 0, 
       0, 0, v*cos(th), 0, sin(th), 0,
       0, 0, 0, (v/wheelBase)*pow(1/cos(psi),2), (1/wheelBase)*tan(psi), 0,
       0, 0, 0, 0, 0, 1,
       0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0;
  Ac = A;
  

  // B matrix
  Eigen::Matrix<double, 6, 3> B;
  B << 0, 0, v*th*sin(th),
       0, 0, -v*th*cos(th),
       0, 0, -(1/wheelBase)*v*psi*pow(1/cos(psi), 2),
       0, 0, 0,
       1, 0, 0,
       0, 1, 0;
  Bc = B;

}

// continuous to discrete matrix conversion
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

int main()
{
    // constants
    const double inf = std::numeric_limits<double>::max();
    const double dT = 0.1;

    /* Ackermann vehicle (LTV) */

    // constants
    const double wheelBase = 1;

    // initial state
    Eigen::Vector<double, 6> x0;
    x0 << 0, 0, 0, 0, 5, 0;

    // linearize
    Eigen::Matrix<double,6,6> Ac;
    Eigen::Matrix<double,6,3> Bc;
    ackermannLinearizedDynMatrices(x0, wheelBase, Ac, Bc);
    std::cout << "Ac = " << Ac << std::endl;
    std::cout << "Bc = " << Bc << std::endl;

    // discretize
    Eigen::MatrixXd Ad;
    Eigen::MatrixXd Bd;
    cont2DiscDynMatrices(Ac, Bc, dT, Ad, Bd);
    std::cout << "Ad = " << Ad << std::endl;
    std::cout << "Bd = " << Bd << std::endl; 

    /* double integrator (LTI) */
    /*

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

    // build MPC object from default constructor
    MpcController MPC;
    MPC.setDynMatrices(A_dyn, B_dyn);
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
      u = MPC.control(x, x_ref);

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
    
    */

    return 0;
}