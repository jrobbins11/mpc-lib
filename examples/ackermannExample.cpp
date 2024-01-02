#include <iostream>
#include <fstream>
#include <limits>
#include <ctime>
#include <cmath>
#include "MpcController.hpp"
#include "MpcMathUtilities.hpp"

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

// function to construct a reference
Eigen::MatrixXd constructReference(const Eigen::VectorXd &x0, int n_horizon, unsigned int n_sim, double dT)
{
  // constants
  const double pi = M_PI;

  // initialize output
  Eigen::MatrixXd x_ref = Eigen::MatrixXd::Zero(6, n_horizon);

  // rotational and velocity frequencies
  double om_th = 2*pi/50;
  double om_v = 2*pi/10;

  // velocity oscillation mean and magnitude
  double v_avg = 5;
  double v_amp = 1;

  // integrate to current sim step
  Eigen::VectorXd x_n = x0; // init
  double t, th, v;
  for (int j=0; j<n_sim; j++)
  {
    // get v and th
    t = j*dT;
    th = x0[2] + om_th*t;
    v = x0[4] + v_amp*sin(om_v*t);

    // integrate
    x_n(0) = x_n(0) + v*cos(th)*dT;
    x_n(1) = x_n(1) + v*sin(th)*dT;
    x_n(2) = th;
    x_n(3) = 0;
    x_n(4) = v;
    x_n(5) = 0;
  }

  // integrate over horizon
  for (int k=0; k<n_horizon; k++)
  {

    // get v and th
    t = (n_sim+k)*dT;
    th = x0[2] + 2*pi*om_th*t;
    v = x0[4] + v_amp*sin(om_v*t);

    // integrate
    if (k==0)
    {
        x_ref(0,k) = x_n(0) + v*cos(th)*dT;
        x_ref(1,k) = x_n(1) + v*sin(th)*dT;
    }
    else
    {
        x_ref(0,k) = x_ref(0,k-1) + v*cos(th)*dT;
        x_ref(1,k) = x_ref(1,k-1) + v*sin(th)*dT;
    }
    x_ref(2,k) = th;
    x_ref(3,k) = 0;
    x_ref(4,k) = v;
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

    // initial state
    Eigen::Vector<double, 6> x0;
    x0 << 0, 0, 0, 0, 5, 0;

    // initial reference
    Eigen::MatrixXd x_ref;
    x_ref = constructReference(x0, n_horizon, 0, dT);

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

    // declare position error variable
    double e_pos;

    // time execution
    clock_t timer;
    float avg_time, time;
    float tot_time = 0;
    float max_time = 0;

    // log sim data
    std::ofstream refFile("../data/ackermannRefData.txt");
    std::ofstream stateFile("../data/ackermannStateData.txt");
    std::ofstream inputFile("../data/ackermannInputData.txt");
    std::ofstream tFile("../data/ackermannTimeData.txt");

    // run MPC controller in loop and print control inputs to screen
    int n_cycles = 500; // number of cycles 
    Eigen::Vector<double, 6> x = x0;
    for (int i=0; i<n_cycles; i++)
    {
     
      // generate reference
      x_ref = constructReference(x0, n_horizon, i, dT);

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

      // step dynamics
      x = ackermannDyn(x, u.segment(0,2), wheelBase, dT);

      // print instantaneous position error magnitude
      // e_pos = sqrt(pow(x(0)-x_ref(0,0), 2) + (pow(x(1)-x_ref(1,0), 2)));
      // std::cout << "e_pos = " << e_pos << std::endl;

      // data logging
      tFile << i*dT << std::endl;;

      for (int j=0; j<6; j++)
        stateFile << x(j) << " ";
      stateFile << std::endl;

      for (int k=0; k<n_horizon; k++)
      {
        for (int j=0; j<6; j++)
          refFile << x_ref(j,k) << " ";
        refFile << std::endl;
      }
      refFile << std::endl;

      for (int j=0; j<2; j++)
        inputFile << u(j) << " ";
      inputFile << std::endl;

    }
    
    // display execution time
    std::cout << "Average MPC execution time = " << 
        tot_time/n_cycles << " sec" << std::endl;
    std::cout << "Max MPC execution time = " <<
        max_time << " sec" << std::endl;

    return 0;
}