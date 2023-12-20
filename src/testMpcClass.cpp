#include <iostream>
#include <limits>
#include <ctime>
#include <cmath>
#include "MpcController.hpp"

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

// function to construct a reference
Eigen::MatrixXd constructReference(const Eigen::VectorXd &x0, double v_des, double th_des, int n_horizon, double t, double dT)
{
  // initialize output
  Eigen::MatrixXd x_ref = Eigen::MatrixXd::Zero(6, n_horizon);

  // construct reference
  for (int k=0; k<n_horizon; k++)
  {
    if (k==0)
    {
      x_ref(0,k) = (x0(0) + v_des*cos(th_des)*t) + v_des*cos(th_des)*dT;
      x_ref(1,k) = (x0(0) + v_des*sin(th_des)*t) + v_des*sin(th_des)*dT;
    }
    else
    {
      x_ref(0,k) = x_ref(0,k-1) + v_des*cos(th_des)*dT;
      x_ref(1,k) = x_ref(1,k-1) + v_des*sin(th_des)*dT;
    }
    x_ref(2,k) = th_des;
    x_ref(3,k) = 0;
    x_ref(4,k) = v_des;
    x_ref(5,k) = 0; 
  }
  return x_ref;
}

int main()
{
    // constants
    const double inf = std::numeric_limits<double>::max();
    const double pi = M_PI;

    // sim parameters
    const double dT = 0.1;
    const int n_horizon = 20;

    /* Ackermann vehicle (LTV) */

    // parameters
    const double wheelBase = 1;

    // state constraints
    Eigen::Matrix<double, 6, 6> I6 = Eigen::Matrix<double, 6, 6>::Identity();
    Eigen::Matrix<double, 6, 6> Ax_ineq;
    Eigen::Vector<double, 6> bx_ineq_low;
    Eigen::Vector<double, 6> bx_ineq_up;
    Ax_ineq << I6;
    bx_ineq_low << -1e4, -1e4, -1e4, -pi/4, -5, -pi;
    bx_ineq_up << 1e4, 1e4, 1e4, pi/4, 10, pi;

    // input constraints
    Eigen::Matrix<double, 3, 3> I3 = Eigen::Matrix<double, 3, 3>::Identity();
    Eigen::Matrix<double, 3, 3> Au_ineq;
    Eigen::Vector<double, 3> bu_ineq_low;
    Eigen::Vector<double, 3> bu_ineq_up;
    Au_ineq << I3;
    bu_ineq_low << -10, -10*pi, 1;
    bu_ineq_up << 10, 10*pi, 1;   

    // cost function matrices
    Eigen::DiagonalMatrix<double, 6> Q;
    Eigen::DiagonalMatrix<double, 3> R;
    Eigen::DiagonalMatrix<double, 6> P;

    Q.diagonal() << 1, 1, 0, 0, 0, 0;
    R.diagonal() << 1e-3, 1e-3*pow(180/pi,2), 0;
    P.diagonal() << 1, 1, 1*pow(180/pi,2), 0, 0, 0;

    // desired velocity and heading
    double v_des, th_des;
    v_des = 7;
    th_des = pi/4;   

    // initial state
    Eigen::Vector<double, 6> x0;
    x0 << 0, 0, 0, 0, 5, 0;

    // initial reference
    Eigen::MatrixXd x_ref;
    x_ref = constructReference(x0, v_des, th_des, n_horizon, 0, dT);

    // initial dynamic linearization (about reference)
    std::vector<Eigen::MatrixXd> A_vec;
    std::vector<Eigen::MatrixXd> B_vec;
      
    Eigen::Matrix<double,6,6> Ac;
    Eigen::Matrix<double,6,3> Bc;
    Eigen::MatrixXd Ad, Bd; // need to be delcared as MatrixXd for cont2DiscDynMatrices

    for (int k=0; k<n_horizon; k++)
    {
      // linearize
      ackermannLinearizedDynMatrices(x_ref.block(0,k,6,1), wheelBase, Ac, Bc);

      // discretize
      cont2DiscDynMatrices(Ac, Bc, dT, Ad, Bd);
      A_vec.push_back(Ad);
      B_vec.push_back(Bd);
    }

    // build MPC object from default constructor
    MpcController MPC;
    MPC.setDynMatrices(A_vec, B_vec);
    MPC.setStageCost(Q, R);
    MPC.setTerminalCost(P);
    MPC.setStateConstraints(Ax_ineq, bx_ineq_low, bx_ineq_up);
    MPC.setInputConstraints(Au_ineq, bu_ineq_low, bu_ineq_up);
    MPC.setMpcHorizon(n_horizon);
    MPC.setMaxExecutionTime(dT);
    MPC.buildController();

    // declare control input variable
    Eigen::Vector<double, 3> u; 

    // time execution
    clock_t timer;
    float avg_time, time;
    float tot_time = 0;
    float max_time = 0;

    // run MPC controller in loop and print control inputs to screen
    int n_cycles = 30; // numbe of cycles 
    Eigen::Vector<double, 6> x = x0;
    for (int i=0; i<n_cycles; i++)
    {
     
      // generate reference
      x_ref = constructReference(x0, v_des, th_des, n_horizon, dT*i, dT);

      // linearize about reference
      for (int k=0; k<n_horizon; k++)
      {
        // linearize
        ackermannLinearizedDynMatrices(x_ref.block(0,k,6,1), wheelBase, Ac, Bc);

        // discretize
        cont2DiscDynMatrices(Ac, Bc, dT, Ad, Bd);
        A_vec[k] = Ad;
        B_vec[k] = Bd;
      }

      // start timer
      timer = clock();

      // compute control
      u = MPC.control(x, x_ref, A_vec, B_vec);

      // log execution time
      timer = clock() - timer;
      time = (float) timer/CLOCKS_PER_SEC;
      tot_time += time;
      if (time > max_time)
        max_time = time;

      // print control output to console
      std::cout << "u = " << u << std::endl;

      // step dynamics
      x = ackermannDyn(x, u.segment(0,2), wheelBase, dT);
    }
    
    // display execution time
    std::cout << "Average MPC execution time = " << 
        tot_time/n_cycles << " sec" << std::endl;
    std::cout << "Max MPC execution time = " <<
        max_time << " sec" << std::endl;



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
      x = A_dyn*x + B_dyn*u;
    }
    
    // display execution time
    std::cout << "Average MPC execution time = " << 
        tot_time/n_cycles << " sec" << std::endl;
    std::cout << "Max MPC execution time = " <<
        max_time << " sec" << std::endl;
    
    */

    return 0;
}