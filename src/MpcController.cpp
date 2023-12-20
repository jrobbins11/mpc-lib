#include "MpcController.hpp"

// problem setup methods
void MpcController::setDynMatrices(const Eigen::MatrixXd &A_dyn_in, 
    const Eigen::MatrixXd &B_dyn_in)
{
    // copy in dynamics matrices
    A_dyn = A_dyn_in;
    B_dyn = B_dyn_in;

    // set state and input dimensions
    n_states = A_dyn.rows();
    n_inputs = B_dyn.cols();
}

void MpcController::setDynMatrices(const std::vector<Eigen::MatrixXd> &A_dyn_vec_in, const std::vector<Eigen::MatrixXd> &B_dyn_vec_in)
{
    // LTV method -> this is an LTV problem
    LTV = true;

    // copy in matrices
    A_dyn_vec = A_dyn_vec_in;
    B_dyn_vec = B_dyn_vec_in;

    // set LTI dynamics matrices (unused)
    A_dyn = A_dyn_vec[0];
    B_dyn = B_dyn_vec[0];

    // set state and input dimensions
    n_states = A_dyn_vec[0].rows();
    n_inputs = B_dyn_vec[0].cols();

}

void MpcController::setStageCost(const Eigen::MatrixXd &Q_cost_in, 
    const Eigen::MatrixXd &R_cost_in)
{
    // copy in cost matrices
    Q_cost = Q_cost_in;
    R_cost = R_cost_in;

    // if terminal cost not already specified, use Q matrix
    if (!terminalCostSpecified)
        P_cost = Q_cost_in;
}

void MpcController::setTerminalCost(const Eigen::MatrixXd &P_cost_in)
{
    // copy in terminal cost matrix
    P_cost = P_cost_in;

    // set terminal cost flag
    terminalCostSpecified = true;
}

void MpcController::setStateConstraints(const Eigen::MatrixXd &Ax_ineq_in, 
      const Eigen::MatrixXd &bx_ineq_low_in, const Eigen::MatrixXd &bx_ineq_up_in)
{
    // copy in matrices
    Ax_ineq = Ax_ineq_in;
    bx_ineq_low = bx_ineq_low_in;
    bx_ineq_up = bx_ineq_up_in;

    // set constraint dimension
    n_xineq = Ax_ineq.rows();

    // if terminal constraint not already specified, use these matrices
    if (!terminalStateConstraintSpecified)
    {
        Ax_term_ineq = Ax_ineq_in;
        bx_term_ineq_low = bx_ineq_low_in;
        bx_term_ineq_up = bx_ineq_up_in;

        // terminal constraint dimension
        n_xtermineq = Ax_term_ineq.rows();
    }

    // constraint flag
    stateConstraintSpecified = true;
}

void MpcController::setInputConstraints(const Eigen::MatrixXd &Au_ineq_in, 
      const Eigen::MatrixXd &bu_ineq_low_in, const Eigen::MatrixXd &bu_ineq_up_in)
{
    // copy in matrices
    Au_ineq = Au_ineq_in;
    bu_ineq_low = bu_ineq_low_in;
    bu_ineq_up = bu_ineq_up_in;

    // input constraint dimension
    n_uineq = Au_ineq.rows();

    // constraint flag
    inputConstraintSpecified = true;
}

void MpcController::setTerminalStateConstraints(const Eigen::MatrixXd &Ax_term_ineq_in, 
      const Eigen::MatrixXd &bx_term_ineq_low_in, const Eigen::MatrixXd &bx_term_ineq_up_in)
{
    // copy in matrices
    Ax_term_ineq = Ax_term_ineq_in;
    bx_term_ineq_low = bx_term_ineq_low_in;
    bx_term_ineq_up = bx_term_ineq_up_in;

    // terminal state constraint dimension
    n_xtermineq = Ax_term_ineq.rows();

    // constraint flag
    terminalStateConstraintSpecified = true;
}

void MpcController::setStateConstraintSlackVarsCost(const Eigen::MatrixXd &Qx_constraint_cost_in)
{
    // copy in matrix
    Qx_constraint_cost = Qx_constraint_cost_in;

    // if terminal constraint softening cost is not specified, use this matrix
    if (!terminalStateConstraintSofteningCostSpecified)
        Qxterm_constraint_cost = Qx_constraint_cost_in;

    // softened state constraints flag
    softenedStateConstraints = true;
}

void MpcController::setInputConstraintSlackVarsCost(const Eigen::MatrixXd &Qu_constraint_cost_in)
{
    // copy in matrix
    Qu_constraint_cost = Qu_constraint_cost_in;

    // softened input constraints flag
    softenedInputConstraints = true;
}

void MpcController::setStateTerminalConstraintSlackVarsCost(const Eigen::MatrixXd &Qxterm_constraint_cost_in)
{
    // copy in matrix
    Qxterm_constraint_cost = Qxterm_constraint_cost_in;

    // set flags
    terminalStateConstraintSofteningCostSpecified = true;
    softenedStateConstraints = true;
}

void MpcController::setMpcHorizon(int n_horizon_in)
{
    n_horizon = n_horizon_in;
}

void MpcController::setMaxExecutionTime(double T_loop_sec_in)
{
    T_sec = T_loop_sec_in;
}

// method to build controller
void MpcController::buildController()
{
    // if constraints aren't specified, build them up as +/- inf
    if (!stateConstraintSpecified)
    {
        Ax_ineq = Eigen::MatrixXd::Zero(1, n_states);
        bx_ineq_low = Eigen::VectorXd::Zero(1);
        bx_ineq_low(0) = -inf;
        bx_ineq_up = Eigen::VectorXd::Zero(1);
        bx_ineq_up(0) = inf;
        n_xineq = 1;

        if (!terminalStateConstraintSpecified)
        {
            Ax_term_ineq = Eigen::MatrixXd::Zero(1, n_states);
            bx_term_ineq_low = Eigen::VectorXd::Zero(1);
            bx_term_ineq_low(0) = -inf;
            bx_term_ineq_up = Eigen::VectorXd::Zero(1);
            bx_term_ineq_up(0) = inf;
            n_xtermineq = 1;
        }
    }

    if (!inputConstraintSpecified)
    {
        Au_ineq = Eigen::MatrixXd::Zero(1, n_inputs);
        bu_ineq_low = Eigen::VectorXd::Zero(1);
        bu_ineq_low(0) = -inf;
        bu_ineq_up = Eigen::VectorXd::Zero(1);
        bu_ineq_up(0) = inf;
        n_uineq = 1;
    }

    // initialize x0 and x_ref
    x0 = Eigen::VectorXd::Zero(n_states);
    x_ref = Eigen::MatrixXd::Zero(n_states, n_horizon);

    // check input dimensions
    validityCheckDimensions();

    // create matrices for optimization problem
    makeInequalityMatrices();
    makeCostMatrices();

    // settings
    solver.settings()->setVerbosity(false);
    solver.settings()->setWarmStart(true);
    solver.settings()->setTimeLimit(T_sec);

    // initial solver setup and configuration
    if (!configureSolver())
        throw std::invalid_argument("Unable to configure solver"); // TO DO: replace with fault behavior?
}

// constructor - no constraint softening
MpcController::MpcController(
    const Eigen::MatrixXd &A_dyn_in, 
    const Eigen::MatrixXd &B_dyn_in,
    const Eigen::MatrixXd &Q_cost_in, 
    const Eigen::MatrixXd &R_cost_in, 
    const Eigen::MatrixXd &P_cost_in, 
    const Eigen::MatrixXd &Ax_ineq_in, 
    const Eigen::VectorXd &bx_ineq_low_in, 
    const Eigen::VectorXd &bx_ineq_up_in,
    const Eigen::MatrixXd &Ax_term_ineq_in, 
    const Eigen::VectorXd &bx_term_ineq_low_in, 
    const Eigen::VectorXd &bx_term_ineq_up_in,
    const Eigen::MatrixXd &Au_ineq_in, 
    const Eigen::VectorXd &bu_ineq_low_in,
    const Eigen::VectorXd &bu_ineq_up_in, 
    int n_horizon_in, 
    double T_loop_sec_in)
{

    // call setup methods
    setDynMatrices(A_dyn_in, B_dyn_in);
    setStageCost(Q_cost_in, R_cost_in);
    setTerminalCost(P_cost_in);
    setStateConstraints(Ax_ineq_in, bx_ineq_low_in, bx_ineq_up_in);
    setTerminalStateConstraints(Ax_term_ineq_in, bx_term_ineq_low_in, bx_term_ineq_up_in);
    setInputConstraints(Au_ineq_in, bu_ineq_low_in, bu_ineq_up_in);
    setMpcHorizon(n_horizon_in);
    setMaxExecutionTime(T_loop_sec_in);

    // build controller
    buildController();

}

// constructor - softened state constraints
MpcController::MpcController(
    const Eigen::MatrixXd &A_dyn_in, 
    const Eigen::MatrixXd &B_dyn_in,
    const Eigen::MatrixXd &Q_cost_in, 
    const Eigen::MatrixXd &R_cost_in, 
    const Eigen::MatrixXd &P_cost_in, 
    const Eigen::MatrixXd &Qx_constraint_cost_in,
    const Eigen::MatrixXd &Qxterm_constraint_cost_in,
    const Eigen::MatrixXd &Ax_ineq_in, 
    const Eigen::VectorXd &bx_ineq_low_in, 
    const Eigen::VectorXd &bx_ineq_up_in, 
    const Eigen::MatrixXd &Ax_term_ineq_in, 
    const Eigen::VectorXd &bx_term_ineq_low_in, 
    const Eigen::VectorXd &bx_term_ineq_up_in,
    const Eigen::MatrixXd &Au_ineq_in, 
    const Eigen::VectorXd &bu_ineq_low_in, 
    const Eigen::VectorXd &bu_ineq_up_in,
    int n_horizon_in, 
    double T_loop_sec_in)
{

    // call setup methods
    setDynMatrices(A_dyn_in, B_dyn_in);
    setStageCost(Q_cost_in, R_cost_in);
    setTerminalCost(P_cost_in);
    setStateConstraintSlackVarsCost(Qx_constraint_cost_in);
    setStateTerminalConstraintSlackVarsCost(Qxterm_constraint_cost_in);
    setStateConstraints(Ax_ineq_in, bx_ineq_low_in, bx_ineq_up_in);
    setTerminalStateConstraints(Ax_term_ineq_in, bx_term_ineq_low_in, bx_term_ineq_up_in);
    setInputConstraints(Au_ineq_in, bu_ineq_low_in, bu_ineq_up_in);
    setMpcHorizon(n_horizon_in);
    setMaxExecutionTime(T_loop_sec_in);

    // build controller
    buildController();


}

// constructor - softened state and input constraints
MpcController::MpcController(
    const Eigen::MatrixXd &A_dyn_in, 
    const Eigen::MatrixXd &B_dyn_in,
    const Eigen::MatrixXd &Q_cost_in, 
    const Eigen::MatrixXd &R_cost_in, 
    const Eigen::MatrixXd &P_cost_in, 
    const Eigen::MatrixXd &Qx_constraint_cost_in, 
    const Eigen::MatrixXd &Qxterm_constraint_cost_in, 
    const Eigen::MatrixXd &Qu_constraint_cost_in,
    const Eigen::MatrixXd &Ax_ineq_in, 
    const Eigen::VectorXd &bx_ineq_low_in, 
    const Eigen::VectorXd &bx_ineq_up_in, 
    const Eigen::MatrixXd &Ax_term_ineq_in, 
    const Eigen::VectorXd &bx_term_ineq_low_in, 
    const Eigen::VectorXd &bx_term_ineq_up_in,
    const Eigen::MatrixXd &Au_ineq_in, 
    const Eigen::VectorXd &bu_ineq_low_in, 
    const Eigen::VectorXd &bu_ineq_up_in,
    int n_horizon_in, 
    double T_loop_sec_in)
{

    // call setup methods
    setDynMatrices(A_dyn_in, B_dyn_in);
    setStageCost(Q_cost_in, R_cost_in);
    setTerminalCost(P_cost_in);
    setStateConstraintSlackVarsCost(Qx_constraint_cost_in);
    setStateTerminalConstraintSlackVarsCost(Qxterm_constraint_cost_in);
    setInputConstraintSlackVarsCost(Qu_constraint_cost_in);
    setStateConstraints(Ax_ineq_in, bx_ineq_low_in, bx_ineq_up_in);
    setTerminalStateConstraints(Ax_term_ineq_in, bx_term_ineq_low_in, bx_term_ineq_up_in);
    setInputConstraints(Au_ineq_in, bu_ineq_low_in, bu_ineq_up_in);
    setMpcHorizon(n_horizon_in);
    setMaxExecutionTime(T_loop_sec_in);

    // build controller
    buildController();

}

// MPC problem setup
void MpcController::makeInequalityMatrices()
{
    // optimization vector z = [x0, ..., xN, u0, ..., uN-1]

    // identity matrices
    Eigen::MatrixXd I_st = Eigen::MatrixXd::Identity(n_states, n_states); // state dimension
    Eigen::MatrixXd I_xineq = Eigen::MatrixXd::Identity(n_xineq, n_xineq); // state constraint dimension
    Eigen::MatrixXd I_xtermineq = Eigen::MatrixXd::Identity(n_xtermineq, n_xtermineq); // state constraint dimension
    Eigen::MatrixXd I_uineq = Eigen::MatrixXd::Identity(n_uineq, n_uineq); // input constraint dimension
    
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
            if (LTV)
            {
                Aeq_states.block(i*n_states,(i-1)*n_states,n_states,n_states) = A_dyn_vec[i-1];
                Aeq_inputs.block(i*n_states,(i-1)*n_inputs,n_states,n_inputs) = B_dyn_vec[i-1];
            }
            else
            {
                Aeq_states.block(i*n_states,(i-1)*n_states,n_states,n_states) = A_dyn;
                Aeq_inputs.block(i*n_states,(i-1)*n_inputs,n_states,n_inputs) = B_dyn;
            }
        }
    }

    // matrix sizes
    int m_Aineq_states = n_horizon*n_xineq;
    int n_Aineq_states = n_horizon*n_states;
    int m_Aineq_states_term = n_xtermineq;
    int n_Aineq_states_term = n_states;
    int m_Aineq_inputs = n_horizon*n_uineq;
    int n_Aineq_inputs = n_horizon*n_inputs;
    int m_b_states = m_Aineq_states;
    int m_b_states_term = m_Aineq_states_term;
    int m_b_inputs = m_Aineq_inputs;

    int m_Aineq_slack_states, n_Aineq_slack_states, m_Aineq_slack_states_term, n_Aineq_slack_states_term;
    if (softenedStateConstraints)
    {
        m_Aineq_slack_states = n_horizon*n_xineq;
        n_Aineq_slack_states = n_horizon*n_xineq;
        m_Aineq_slack_states_term = n_xtermineq;
        n_Aineq_slack_states_term = n_xtermineq;
    }
    else
    {
        m_Aineq_slack_states = 0;
        n_Aineq_slack_states = 0;
        m_Aineq_slack_states_term = 0;
        n_Aineq_slack_states_term = 0;
    }

    int m_Aineq_slack_inputs, n_Aineq_slack_inputs;
    if (softenedInputConstraints)
    {
        m_Aineq_slack_inputs = n_horizon*n_uineq;
        n_Aineq_slack_inputs = n_horizon*n_uineq;
    }
    else
    {
        m_Aineq_slack_inputs = 0;
        n_Aineq_slack_inputs = 0;
    }
    
    // setup state and input inequality constraints, initialize to zero
    Eigen::MatrixXd Aineq_states = Eigen::MatrixXd::Zero(m_Aineq_states, n_Aineq_states);
    Eigen::MatrixXd Aineq_states_term = Eigen::MatrixXd::Zero(m_Aineq_states_term, n_Aineq_states_term);
    Eigen::MatrixXd Aineq_inputs = Eigen::MatrixXd::Zero(m_Aineq_inputs, n_Aineq_inputs);
    Eigen::MatrixXd Aineq_slack_states = Eigen::MatrixXd::Zero(m_Aineq_slack_states, n_Aineq_slack_states);
    Eigen::MatrixXd Aineq_slack_states_term = Eigen::MatrixXd::Zero(m_Aineq_slack_states_term, n_Aineq_slack_states_term);
    Eigen::MatrixXd Aineq_slack_inputs = Eigen::MatrixXd::Zero(m_Aineq_slack_inputs, n_Aineq_slack_inputs);
    Eigen::VectorXd bup_states = Eigen::VectorXd::Zero(m_b_states);
    Eigen::VectorXd blow_states = Eigen::VectorXd::Zero(m_b_states);
    Eigen::VectorXd bup_states_term = Eigen::VectorXd::Zero(m_b_states_term);
    Eigen::VectorXd blow_states_term = Eigen::VectorXd::Zero(m_b_states_term);
    Eigen::VectorXd bup_inputs = Eigen::VectorXd::Zero(m_b_inputs);
    Eigen::VectorXd blow_inputs = Eigen::VectorXd::Zero(m_b_inputs);

    // write state inequalities
    for (int i=0; i<n_horizon; i++) 
    {   
        
        if (i==0) // no constraint imposed on x0
        {
            for (int j=0; j<n_xineq; j++)
            {
                blow_states(i*n_xineq+j) = -inf;
                bup_states(i*n_xineq+j) = inf;
            }
        }
        else 
        {
            Aineq_states.block(i*n_xineq,i*n_states,n_xineq,n_states) = Ax_ineq;
            blow_states.segment(i*n_xineq,n_xineq) = bx_ineq_low;
            bup_states.segment(i*n_xineq,n_xineq) = bx_ineq_up;

            if (softenedStateConstraints)
                Aineq_slack_states.block(i*n_xineq,i*n_xineq,n_xineq,n_xineq) = -I_xineq;
        }
    }

    // write terminal state inequalities
    Aineq_states_term = Ax_term_ineq;
    blow_states_term = bx_term_ineq_low;
    bup_states_term = bx_term_ineq_up;
    if (softenedStateConstraints)
        Aineq_slack_states_term = -I_xtermineq;

    // write input inequalities
    for (int i=0; i<n_horizon; i++)
    {
        Aineq_inputs.block(i*n_uineq, i*n_inputs, n_uineq, n_inputs) = Au_ineq;
        blow_inputs.segment(i*n_uineq, n_uineq) = bu_ineq_low;
        bup_inputs.segment(i*n_uineq, n_uineq) = bu_ineq_up;

        if (softenedInputConstraints)
            Aineq_slack_inputs.block(i*n_uineq,i*n_uineq,n_uineq,n_uineq) = -I_uineq;
    }

    // resize (and init to zero) inequality matrices
    A.resize(m_Aeq_states + m_Aineq_states + m_Aineq_states_term + m_Aineq_inputs, 
        n_Aeq_states + n_Aeq_inputs + n_Aineq_slack_states + n_Aineq_slack_states_term + n_Aineq_slack_inputs);
    b_low = Eigen::VectorXd::Zero(m_beq + m_b_states + m_b_states_term + m_b_inputs);
    b_up = Eigen::VectorXd::Zero(m_beq + m_b_states + m_b_states_term + m_b_inputs);

    // get triplets to fill sparse matrix
    std::vector<Eigen::Triplet<double>> tripvec;
    getTripletsForMatrix(Aeq_states, tripvec, 0, 0);
    getTripletsForMatrix(Aeq_inputs, tripvec, 0, n_Aeq_states);
    getTripletsForMatrix(Aineq_states, tripvec, m_Aeq_states, 0);
    getTripletsForMatrix(Aineq_states_term, tripvec, m_Aeq_states+m_Aineq_states, n_Aineq_states);
    getTripletsForMatrix(Aineq_inputs, tripvec, m_Aeq_states+m_Aineq_states+m_Aineq_states_term, 
        n_Aineq_states+n_Aineq_states_term);
    if (softenedStateConstraints)
        getTripletsForMatrix(Aineq_slack_states, tripvec, m_Aeq_states, n_Aeq_states+n_Aeq_inputs);
        getTripletsForMatrix(Aineq_slack_states_term, tripvec, m_Aeq_states+m_Aineq_states, 
            n_Aeq_states+n_Aeq_inputs+n_Aineq_slack_states);
    if (softenedInputConstraints)
        getTripletsForMatrix(Aineq_slack_inputs, tripvec, m_Aeq_states+m_Aineq_states+m_Aineq_states_term, 
            n_Aeq_states+n_Aeq_inputs+n_Aineq_slack_states+n_Aineq_slack_states_term);

    // fill sparse matrix
    A.setFromTriplets(tripvec.begin(), tripvec.end());

    // set b_low and b_up
    b_low.segment(0,m_beq) = beq;
    b_up.segment(0,m_beq) = beq;

    b_low.segment(m_beq,m_b_states) = blow_states;
    b_up.segment(m_beq,m_b_states) = bup_states;

    b_low.segment(m_beq+m_b_states,m_b_states_term) = blow_states_term;
    b_up.segment(m_beq+m_b_states,m_b_states_term) = bup_states_term;

    b_low.segment(m_beq+m_b_states+m_b_states_term,m_b_inputs) = blow_inputs;
    b_up.segment(m_beq+m_b_states+m_b_states_term,m_b_inputs) = bup_inputs;

}

void MpcController::makeCostMatrices()
{
    // dimensions
    int m_H_states = n_states*(n_horizon+1);
    int n_H_states = m_H_states;
    int m_H_inputs = n_inputs*n_horizon;
    int n_H_inputs = m_H_inputs;

    int m_H_slack_states, n_H_slack_states, m_H_slack_states_term, n_H_slack_states_term;
    if (softenedStateConstraints)
    {
        m_H_slack_states = n_xineq*n_horizon;
        n_H_slack_states = m_H_slack_states;
        m_H_slack_states_term = n_xtermineq;
        n_H_slack_states_term = m_H_slack_states_term;
    }
    else
    {
        m_H_slack_states = 0;
        n_H_slack_states = 0;
        m_H_slack_states_term = 0;
        n_H_slack_states_term = 0;
    }

    int m_H_slack_inputs, n_H_slack_inputs;
    if (softenedInputConstraints)
    {
        m_H_slack_inputs = n_uineq*n_horizon;
        n_H_slack_inputs = m_H_slack_inputs;
    }
    else
    {
        m_H_slack_inputs = 0;
        n_H_slack_inputs = 0;
    }
    
    // resize matrices
    H.resize(m_H_states + m_H_inputs + m_H_slack_states + m_H_slack_states_term + m_H_slack_inputs, 
        n_H_states + n_H_inputs + n_H_slack_states + n_H_slack_states_term + n_H_slack_inputs);
    f.resize(n_H_states + n_H_inputs + n_H_slack_states + n_H_slack_states_term + n_H_slack_inputs);
    
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
    if (softenedStateConstraints)
    {
        for (int i=1; i<n_horizon; i++) // no cost on x0
        {
            getTripletsForMatrix(Qx_constraint_cost, tripvec, 
                m_H_states + m_H_inputs + n_xineq*i, 
                n_H_states + n_H_inputs + n_xineq*i);
        }
        getTripletsForMatrix(Qxterm_constraint_cost, tripvec,
            m_H_states + m_H_inputs + m_H_slack_states,
            n_H_states + n_H_inputs + n_H_slack_states);
    }
    if (softenedInputConstraints)
    {
        for (int i=0; i<n_horizon; i++)
        {
            getTripletsForMatrix(Qu_constraint_cost, tripvec,
                m_H_states + m_H_inputs + m_H_slack_states + m_H_slack_states_term + n_uineq*i,
                n_H_states + n_H_inputs + n_H_slack_states + n_H_slack_states_term + n_uineq*i);
        }
    }

    // fill H
    H.setFromTriplets(tripvec.begin(), tripvec.end());

    // fill f
    f = Eigen::VectorXd::Zero(n_H_states+n_H_inputs+n_H_slack_states+n_H_slack_states_term+n_H_slack_inputs); // init to zero
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

    // check argument dimensions
    if (x.rows() != n_states)
        throw std::invalid_argument("Inconsistent state vector dimension"); // TO DO: fault behavior

    if ((x_ref.rows() != n_states) || (x_ref.cols() != n_horizon))
        throw std::invalid_argument("Inconsistent state reference vector dimensions"); // TO DO: fault behavior

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

// LTV control method
Eigen::VectorXd MpcController::control(const Eigen::VectorXd &x, const Eigen::MatrixXd &x_ref,
      const std::vector<Eigen::MatrixXd> &A_dyn_vec, const std::vector<Eigen::MatrixXd> &B_dyn_vec)
{

    // set LTV flag
    LTV = true;

    // check argument dimensions
    if (x.rows() != n_states)
        throw std::invalid_argument("Inconsistent state vector dimension"); // TO DO: fault behavior

    if ((x_ref.rows() != n_states) || (x_ref.cols() != n_horizon))
        throw std::invalid_argument("Inconsistent state reference vector dimensions"); // TO DO: fault behavior

    for (int i=0; i<n_horizon; i++)
    {
        if ((A_dyn_vec[i].rows() != n_states) || (A_dyn_vec[i].cols() != n_states))
            throw std::invalid_argument("LTV A matrix dimensions are incosistent");
        if ((B_dyn_vec[i].rows() != n_states) || (B_dyn_vec[i].cols() != n_inputs))
            throw std::invalid_argument("LTV B matrix dimensions are incosistent");
    }
    // TO DO: check length of dynamic matrix vectors

    // store inputs
    this->x0 = x;
    this->x_ref = x_ref;
    this->A_dyn_vec = A_dyn_vec;
    this->B_dyn_vec = B_dyn_vec;

    // constraint update
    // TO DO: split up this function so that only the dynamic portions of matrices need to be updated
    makeInequalityMatrices();

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

void MpcController::validityCheckDimensions()
{
    // check state dimension
    if ((A_dyn.rows() != n_states) || (A_dyn.cols() != n_states) || 
        (B_dyn.rows() != n_states) || (Q_cost.rows() != n_states) || 
        (Q_cost.cols() != n_states) || (P_cost.rows() != n_states) || 
        (P_cost.cols() != n_states) || (Ax_ineq.cols() != n_states) || 
        (Ax_term_ineq.cols() != n_states) || 
        (x0.rows() != n_states) || (x_ref.rows() != n_states))
    {
        throw std::invalid_argument("Inconsistent state dimensions");
    }
    
    // check input dimension
    if ((R_cost.rows() != n_inputs) || (R_cost.cols() != n_inputs) ||
        (Au_ineq.cols() != n_inputs))
    {
       throw std::invalid_argument("Inconsistent input dimensions");
    }

    // check state constraints
    if ((Ax_ineq.rows() != n_xineq) || (bx_ineq_low.rows() != n_xineq) || (bx_ineq_up.rows() != n_xineq))
        throw std::invalid_argument("Inconsistent state constraint dimensions");
    
    // check terminal state constraints
    if ((Ax_term_ineq.rows() != n_xtermineq) || (bx_term_ineq_low.rows() != n_xtermineq) || (bx_term_ineq_up.rows() != n_xtermineq))
        throw std::invalid_argument("Inconsistent terminal state constraint dimensions");

    // check input constraints
    if ((Au_ineq.rows() != n_uineq) || (bu_ineq_low.rows() != n_uineq) || (bu_ineq_up.rows() != n_uineq))
        throw std::invalid_argument("Inconsistent input constraint dimensions");

    // check reference length
    if (x_ref.cols() != n_horizon)
        throw std::invalid_argument("Inconsistent prediction/control horizon");

    if (softenedStateConstraints)
    {
        // check softened state constraint dimension
        if ((Qx_constraint_cost.rows() != n_xineq) || (Qx_constraint_cost.cols() != n_xineq))
            throw std::invalid_argument("Inconsistent state slack variable dimensions");

        // check softened terminal state constraint dimension
        if ((Qxterm_constraint_cost.rows() != n_xtermineq) || (Qxterm_constraint_cost.cols() != n_xtermineq))
            throw std::invalid_argument("Inconsistent terminal state slack variable dimensions");
    }

    if (softenedInputConstraints)
    {
        // check softened input constraint dimension
        if ((Qu_constraint_cost.rows() != n_uineq) || (Qu_constraint_cost.cols() != n_uineq))
            throw std::invalid_argument("Inconsistent input slack variable dimensions");
    }

    // check consistency of matrix dimensions for LTV problems
    // TO DO: check length of dynamic matrix vectors
    if (LTV)
    {
        for (int i=0; i<n_horizon; i++)
        {
            if ((A_dyn_vec[i].rows() != n_states) || (A_dyn_vec[i].cols() != n_states))
                throw std::invalid_argument("LTV A matrix dimensions are incosistent");
            if ((B_dyn_vec[i].rows() != n_states) || (B_dyn_vec[i].cols() != n_inputs))
                throw std::invalid_argument("LTV B matrix dimensions are incosistent");
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

    // dimensions
    std::cout << "Problem dimensions:" << std::endl;
    std::cout << " number of unknowns = " << A.cols() << std::endl;
    std::cout << " number of constraints = " << A.rows() << std::endl;
}

void MpcController::printSolution()
{
    // header
    std::cout << "Optimal solution vector z is:" << std::endl;

    // print z
    std::cout << solution << std::endl;
}
