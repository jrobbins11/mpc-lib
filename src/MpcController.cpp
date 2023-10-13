#include "MpcController.hpp"

// constructor
MpcController::MpcController(
      const Eigen::VectorXd &x0, const Eigen::MatrixXd &x_ref,
      const Eigen::MatrixXd &A_dyn, const Eigen::MatrixXd &B_dyn,
      const Eigen::MatrixXd &Q_cost, const Eigen::MatrixXd &R_cost, 
      const Eigen::MatrixXd &P_cost, const Eigen::MatrixXd &Ax_ineq, 
      const Eigen::VectorXd &bx_ineq, const Eigen::MatrixXd &Au_ineq, 
      const Eigen::VectorXd &bu_ineq, int n_horizon): 
      x0(x0), x_ref(x_ref), A_dyn(A_dyn), B_dyn(B_dyn), Q_cost(Q_cost), 
      R_cost(R_cost), P_cost(P_cost), Ax_ineq(Ax_ineq), bx_ineq(bx_ineq), 
      Au_ineq(Au_ineq), bu_ineq(bu_ineq), n_horizon(n_horizon)
{
    // validity checking dimensions
    int n_work; // working variable
    n_work = A_dyn.rows(); // state dimension
    if ((A_dyn.cols() != n_work) || (B_dyn.rows() != n_work) || 
        (Q_cost.rows() != n_work) || (Q_cost.cols() != n_work) ||
        (P_cost.rows() != n_work) || (P_cost.cols() != n_work) ||
        (Ax_ineq.cols() != n_work) || (x0.rows() != n_work) ||
        (x_ref.rows() != n_work))
    {
        throw std::invalid_argument("Inconsistent state dimensions");
    }
    
    n_work = B_dyn.cols(); // input dimension
    if ((R_cost.rows() != n_work) || (R_cost.cols() != n_work) ||
        (Au_ineq.cols() != n_work))
    {
       throw std::invalid_argument("Inconsistent input dimensions");
    }

    n_work = Ax_ineq.rows(); // state constraints
    if (bx_ineq.rows() != n_work)
        throw std::invalid_argument("Inconsistent state constraint dimensions");

    n_work = Au_ineq.rows(); // input constraints
    if (bu_ineq.rows() != n_work)
        throw std::invalid_argument("Inconsistent input constraint dimensions");

    if (x_ref.cols() != n_horizon) // reference length
        throw std::invalid_argument("Inconsistent prediction/control horizon");

    // get number of states, constraints, and inputs
    n_states = A_dyn.rows();
    n_inputs = B_dyn.cols();
    n_xineq = Ax_ineq.rows();
    n_uineq = Au_ineq.rows();

    // create matrices for optimization problem
    makeInequalityMatrices();
    makeCostMatrices();

    // solver options
    solver.settings()->setVerbosity(false);
    solver.settings()->setWarmStart(true);

    // initial solver setup and configuration
    if (!configureSolver())
        throw std::invalid_argument("Unable to configure solver"); // TO DO: replace with fault behavior?
    if (!solveOptimizationProblem())
        throw std::invalid_argument("Unable to solver optimization problem"); // TO DO: replace with fault behavior?

}

// MPC problem setup
void MpcController::makeInequalityMatrices()
{
    // optimization vector z = [x0, ..., xN, u0, ..., uN-1]

    // identity matrix with dimensions of state
    Eigen::MatrixXd I_st = Eigen::MatrixXd::Identity(n_states, n_states);
    
    // matrix sizes
    int m_Aeq_states = (n_horizon+1)*n_states;
    int n_Aeq_states = m_Aeq_states;
    int m_Aeq_inputs = m_Aeq_states;
    int n_Aeq_inputs = n_horizon*n_inputs;
    int m_beq = m_Aeq_states;

    // setup dynamic constraints: xk+1 = A*xk + B*uk, x0 = x0
    // initialize matrices to zeros
    Eigen::MatrixXd Aeq_states = Eigen::MatrixXd::Zero(m_Aeq_states, n_Aeq_states);
    Eigen::MatrixXd Aeq_inputs = Eigen::MatrixXd::Zero(m_Aeq_inputs, n_Aeq_inputs);
    Eigen::VectorXd beq = Eigen::VectorXd::Zero(m_beq);

    // write Aeq_states, Aeq_inputs, beq
    for (int i=0; i<(n_horizon+1); i++)
    {
        if (i==0) // x0 = x0 constraint
        {
            Aeq_states.block(i*n_states,i*n_states,n_states,n_states) = -1*I_st;
            beq.segment(i*n_states, n_states) = -1*x0;
        }
        else // dynamic constraints
        {
            Aeq_states.block(i*n_states,i*n_states,n_states,n_states) = -1*I_st;
            Aeq_states.block(i*n_states,(i-1)*n_states,n_states,n_states) = A_dyn;
            Aeq_inputs.block(i*n_states,(i-1)*n_inputs,n_states,n_inputs) = B_dyn;
        }
    }

    // matrix sizes
    int m_Aineq_states = (n_horizon+1)*n_xineq;
    int n_Aineq_states = (n_horizon+1)*n_states;
    int m_Aineq_inputs = n_horizon*n_uineq;
    int n_Aineq_inputs = n_horizon*n_inputs;
    int m_b_states = m_Aineq_states;
    int m_b_inputs = m_Aineq_inputs;
    
    // setup state and input inequality constraints, initialize to zero
    Eigen::MatrixXd Aineq_states = Eigen::MatrixXd::Zero(m_Aineq_states, n_Aineq_states);
    Eigen::MatrixXd Aineq_inputs = Eigen::MatrixXd::Zero(m_Aineq_inputs, n_Aineq_inputs);
    Eigen::VectorXd bup_states = Eigen::VectorXd::Zero(m_b_states);
    Eigen::VectorXd blow_states = Eigen::VectorXd::Zero(m_b_states);
    Eigen::VectorXd bup_inputs = Eigen::VectorXd::Zero(m_b_inputs);
    Eigen::VectorXd blow_inputs = Eigen::VectorXd::Zero(m_b_inputs);

    // write state inequalities
    for (int i=0; i<(n_horizon+1); i++) 
    {   
        for (int j=0; j<n_xineq; j++)
            blow_states(i*n_xineq+j) = -1*OsqpEigen::INFTY;
        
        if (i==0) // no constraint imposed on x0
        {
            for (int j=0; j<n_xineq; j++)
                bup_states(i*n_xineq+j) = OsqpEigen::INFTY;
        }
        else 
        {
            Aineq_states.block(i*n_xineq,i*n_states,n_xineq,n_states) = Ax_ineq;
            bup_states.segment(i*n_xineq,n_xineq) = bx_ineq;
        }
    }

    // write input inequalities
    for (int i=0; i<n_horizon; i++)
    {
        for (int j=0; j<n_uineq; j++)
            blow_inputs(i*n_uineq+j) = -1*OsqpEigen::INFTY;
        Aineq_inputs.block(i*n_uineq, i*n_inputs, n_uineq, n_inputs) = Au_ineq;
        bup_inputs.segment(i*n_uineq, n_uineq) = bu_ineq;
    }

    // resize inequality matrices
    A.resize(m_Aeq_states + m_Aineq_states + m_Aineq_inputs, 
        n_Aeq_states + n_Aeq_inputs);
    b_low.resize(m_Aeq_states + m_Aineq_states + m_Aineq_inputs);
    b_up.resize(m_Aeq_states + m_Aineq_states + m_Aineq_inputs);

    // get triplets to fill sparse matrix
    std::vector<Eigen::Triplet<double>> tripvec;
    getTripletsForMatrix(Aeq_states, tripvec, 0, 0);
    getTripletsForMatrix(Aeq_inputs, tripvec, 0, n_Aeq_states);
    getTripletsForMatrix(Aineq_states, tripvec, m_Aeq_states, 0);
    getTripletsForMatrix(Aineq_inputs, tripvec, m_Aeq_states+m_Aineq_states, n_Aeq_states);

    // fill sparse matrix
    A.setFromTriplets(tripvec.begin(), tripvec.end());

    // set b_low and b_up
    b_low.segment(0,m_beq) = beq;
    b_up.segment(0,m_beq) = beq;

    b_low.segment(m_beq,m_b_states) = blow_states;
    b_up.segment(m_beq,m_b_states) = bup_states;

    b_low.segment(m_beq+m_b_states,m_b_inputs) = blow_inputs;
    b_up.segment(m_beq+m_b_states,m_b_inputs) = bup_inputs;

}

void MpcController::makeCostMatrices()
{
    // dimensions
    int m_H_states = n_states*(n_horizon+1);
    int n_H_states = m_H_states;
    int m_H_inputs = n_inputs*n_horizon;
    int n_H_inputs = m_H_inputs;
    
    // resize matrices
    H.resize(m_H_states+n_H_inputs, m_H_states+n_H_inputs);
    f.resize(n_H_states+n_H_inputs);

    // get triplets to fill H
    std::vector<Eigen::Triplet<double>> tripvec;
    for (int i=1; i<n_horizon+1; i++) // no cost on x0
    {
        if (i < n_horizon)
            getTripletsForMatrix(Q_cost, tripvec, n_states*i, n_states*i);
        else
            getTripletsForMatrix(P_cost, tripvec, n_states*i, n_states*i);
    }
    for (int i=0; i<n_horizon; i++)
    {
        getTripletsForMatrix(R_cost, tripvec, m_H_states + n_inputs*i, n_H_states + n_inputs*i);
    }

    // fill H
    H.setFromTriplets(tripvec.begin(), tripvec.end());

    // fill f
    Eigen::VectorXd x_ref_i; // declare
    for (int i=1; i<n_horizon+1; i++)
    {
        // reference at step i, offset b/c x0 included in state
        x_ref_i = x_ref.block(0, i-1, n_states, 1);

        // insert -Q*x_ref into gradient vector
        if (i < n_horizon)
            f.segment(i*n_states, n_states) = -Q_cost*x_ref_i;
        else
            f.segment(i*n_states, n_states) = -P_cost*x_ref_i;
    }

    // fill remainder with zeros
    f.segment(n_H_states, n_H_inputs) = Eigen::VectorXd::Zero(n_H_inputs);
}

// solver functions
bool MpcController::configureSolver()
{
    // pass matrices and vars to solver
    solver.data()->setNumberOfVariables(H.rows());
    solver.data()->setNumberOfConstraints(A.rows());
    if (!solver.data()->setHessianMatrix(H))
        return false;
    if (!solver.data()->setGradient(f))
        return false;
    if (!solver.data()->setLinearConstraintsMatrix(A))
        return false;
    if (!solver.data()->setLowerBound(b_low))
        return false;
    if (!solver.data()->setUpperBound(b_up))
        return false;

    // instantiate the solver
    if (!solver.initSolver())
        return false;

    // return true if successful
    return true;
}

bool MpcController::solveOptimizationProblem()
{
    // solve
    if (solver.solveProblem() != OsqpEigen::ErrorExitFlag::NoError)
        return false;
    solution = solver.getSolution();

    // return true if successful
    return true;
}

// control function
Eigen::VectorXd MpcController::control(const Eigen::VectorXd &x, const Eigen::MatrixXd &x_ref)
{
    // TO DO: call updateGradient, updateBounds
    // solve optimization problem and return solution
    // https://robotology.github.io/osqp-eigen/class_osqp_eigen_1_1_solver.html

    // store x0 and x_ref
    this->x0 = x;
    this->x_ref = x_ref;
    
    // constraint update
    b_low.segment(0, n_states) = -x;
    b_up.segment(0, n_states) = -x;

    // gradient update
    Eigen::VectorXd x_ref_i; // declare
    for (int i=1; i<n_horizon+1; i++)
    {
        // reference at step i, offset b/c x0 included in state
        x_ref_i = x_ref.block(0, i-1, n_states, 1);

        // insert -Q*x_ref into gradient vector
        if (i < n_horizon)
            f.segment(i*n_states, n_states) = -Q_cost*x_ref_i;
        else
            f.segment(i*n_states, n_states) = -P_cost*x_ref_i;
    }

    // call OsqpEigen update methods
    if (!solver.updateBounds(b_low, b_up))
        std::cout << "failed to update bounds" << std::endl; // TO DO: fault behavior
    if (!solver.updateGradient(f))
        std::cout << "failed to update gradient" << std::endl; // TO DO: fault behavior

    // solve the optimization problem
    if (!solveOptimizationProblem())
        std::cout << "solver returned with error" << std::endl; // TO DO: fault behavior
    
    // return the control input for step 0
    return solution.segment(n_states*(n_horizon+1), n_inputs);

}

// utilities
void MpcController::getTripletsForMatrix(const Eigen::MatrixXd &mat, std::vector<Eigen::Triplet<double>> &tripvec,
      int rowOffset, int colOffset)
{
    int m = mat.rows();
    int n = mat.cols();
    for (int i=0; i<m; i++)
    {
        for (int j=0; j<n; j++)
        {
            if (mat(i,j) != 0)
                tripvec.push_back(Eigen::Triplet<double>(i+rowOffset, j+colOffset, mat(i,j)));
        }
    }
}

// print methods
void MpcController::printOptimizationProblem()
{
    // header
    std::cout << "Quadratic programming problem of form:" << std::endl;
    std::cout << " min(z) 0.5*z'*H*z + f*z" << std::endl;
    std::cout << " subject to:" << std::endl;
    std::cout << " b_low <= A*z <= b_up\n" << std::endl;
    std::cout << "Matrices are:\n" << std::endl;

    // cost
    std::cout << "H = " << std::endl;
    std::cout << Eigen::MatrixXd(H) << std::endl;
    std::cout << "f = " << std::endl;
    std::cout << f << std::endl;

    // constraints
    std::cout << "A = " << std::endl;
    std::cout << Eigen::MatrixXd(A) << std::endl;
    std::cout << "b_low = " << std::endl;
    std::cout << b_low << std::endl;
    std::cout << "b_up = " << std::endl;
    std::cout << b_up << std::endl;
}