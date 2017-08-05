#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// TODO: Set the timestep length and duration
// N = 10, dt = .1 --> Worked Fantastic!
// N = 20, dt = .1 --> A bit off, big turns confuse it
// N = 25, dt = .05 --> Barely left the station
size_t N = 10;
double dt = .1;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG (center of gravity) that has a similar radius.
const double Lf = 2.67;

double ref_cte = 0;
double ref_epsi = 0;
double ref_v = 100;

// Pass this as a 1-dimensional vector
int x_start = 0;
int y_start = x_start + N;
int psi_start = y_start + N;
int v_start = psi_start + N;
int cte_start = v_start + N;
int epsi_start = cte_start + N;
int delta_start = epsi_start + N;
int a_start = delta_start + N - 1;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.
    fg[0] = 0; // cost function

    // Reference State Cost
    // TODO: Define the cost related to the reference state
    // and anything you thiink may be beneficial
    for (int i = 0; i < N; i++) {
      // Coeffiencts have a lot of power on how much intention we
      // want the cost function to pay to certain attributes.
      // The higher the value, the lower we want the value.
      // (i.e. lots of attn to cte cause we want low value)
      // fg[0] += 2000 * CppAD::pow(vars[cte_start + i] - ref_cte, 2);
      // fg[0] += 2000 * CppAD::pow(vars[epsi_start + i] - ref_epsi, 2);
      // fg[0] += CppAD::pow(vars[v_start + i] - ref_v, 2);
      AD<double> cte = vars[cte_start + i];
      AD<double> epsi = vars[epsi_start + i];
      AD<double> v = vars[v_start + i];

      fg[0] += 2000 * CppAD::pow(cte, 2) + 2000 * CppAD::pow(epsi, 2) + CppAD::pow(v - ref_v, 2);
    }

    for (int i = 0; i < N - 1; i++) {
      fg[0] += 5 * CppAD::pow(vars[delta_start + i], 2);
      fg[0] += 5 * CppAD::pow(vars[a_start + i], 2);
    }

    for (int i = 0; i < N - 2; i++) {
      fg[0] += 200 * CppAD::pow(vars[delta_start + i + 1] - vars[delta_start + i], 2);
      fg[0] += 10 * CppAD::pow(vars[a_start + i + 1] - vars[a_start + 1], 2);
    }

    //
    // Setup Constraints
    //
    // Note: In this section, you'll setup the model constraints.

    // Initial Constraints
    //
    // We add 1 to each of the starting indices due to cost being
    // located at index 0 of `fg`.
    // This bumps up the position of all the other values.
    fg[1 + x_start] = vars[x_start];
    fg[1 + y_start] = vars[y_start];
    fg[1 + psi_start] = vars[psi_start];
    fg[1 + v_start] = vars[v_start];
    fg[1 + cte_start] = vars[cte_start];
    fg[1 + epsi_start] = vars[epsi_start];

    // Kinematic model
    // The rest of the constraints
    for (int i = 1; i < N; i++) {
      // state/model @ t + 1
      AD<double> x1 = vars[x_start + i];
      AD<double> y1 = vars[y_start + i];
      AD<double> psi1 = vars[psi_start + i];
      AD<double> v1 = vars[v_start + i];
      AD<double> cte1 = vars[cte_start + i];
      AD<double> epsi1 = vars[epsi_start + i];

      // state/model @ t
      AD<double> x0 = vars[x_start + i - 1];
      AD<double> y0 = vars[y_start + i - 1];
      AD<double> psi0 = vars[psi_start + i - 1];
      AD<double> v0 = vars[v_start + i - 1];
      AD<double> cte0 = vars[cte_start + i - 1];
      AD<double> epsi0 = vars[epsi_start + i - 1];

      AD<double> delta0 = vars[delta_start + i - 1];
      AD<double> a0 = vars[a_start + i - 1];

      // calulate polynomial value w/ coeffs
      AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * CppAD::pow(x0,2) + coeffs[3] * CppAD::pow(x0,3);
      // desired psi value
      AD<double> psides0 = CppAD::atan(coeffs[1] + 2 * coeffs[2] * x0 + 3 * coeffs[3] * CppAD::pow(x0,2));

      // Here's `x` to get you started.
      // The idea here is to constraint this value to be 0.
      //
      // NOTE: The use of `AD<double>` and use of `CppAD`!
      // This is also CppAD can compute derivatives and pass
      // these to the solver.

      // TODO: Setup the rest of the model constraints
      // Model equations shown in lesson/classroom
      // x1 = x @ t + 1
      // x0 = x @ t
      // v0 * CppAD::cos(psi0) * dt = how x0 is changing/rate of change
      fg[1 + x_start + i] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[1 + y_start + i] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[1 + psi_start + i] = psi1 - (psi0 + v0 / Lf * dt * delta0);
      fg[1 + v_start + i] = v1 - (v0 + a0 * dt);
      fg[1 + cte_start + i] = cte1 - (f0 - y0 + v0 * CppAD::sin(epsi0) * dt);
      fg[1 + epsi_start + i] = epsi1 - (psi0 - psides0 + v0 / Lf * delta0 * dt);
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  typedef CPPAD_TESTVECTOR(double) Dvector;
  double x = state[0];
  double y = state[1];
  double psi = state[2];
  double v = state[3];
  double cte = state[4];
  double epsi = state[5];

  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // # number of timesteps
  // 6 => x, y, psi, v. cte, epsi
  // (N - 1) * 2 => delta & acceleration
  // (N - 1) => You don't need delta at the end.
  // Don't know where your're going after
  // size_t n_vars = N * 6 + (N - 1) * 2;
  size_t n_vars = 6 * N + 2 * (N - 1);
  // TODO: Set the number of constraints
  size_t n_constraints = N * 6; // 60 n_constraints because N is 10.

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // TODO: Set lower and upper limits for variables.

  // Set all non-actuators upper and lowerlimits
  // to the max negative and positive values.
  for (int i = 0; i < delta_start; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  // Acceleration/decceleration upper and lower limits.
  // NOTE: Feel free to change this to something else.
  for (int i = a_start; i < n_vars; i++) {
    vars_lowerbound[i] = -1;
    vars_upperbound[i] = 1;
  }

  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  // NOTE: Feel free to change this to something else.
  for (int i = delta_start; i < a_start; i++) {
    vars_lowerbound[i] = -1;
    vars_upperbound[i] = 1;
  }

  // Lower and upper limits for constraints
  // All of these should be 0 except the initial
  // state indices.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  // We make a special case for our first state.
  // So we set all those values to the inital state.
  // If we don't do that, Solver doesnt really know where to start.
  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[psi_start] = psi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[epsi_start] = epsi;

  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[psi_start] = psi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[epsi_start] = epsi;

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  false        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  vector<double> result;


  // Steering and Acutation needed for next time step.
  result.push_back(solution.x[delta_start]);
  result.push_back(solution.x[a_start]);

  for (int i = 0; i < N - 1; i++) {
    // Capture where the MPC will be in the future.
    // Grabbing solutions x & y coordinates for where mpc
    // controller is going. That is the green line.
    result.push_back(solution.x[x_start + i + 1]);
    result.push_back(solution.x[y_start + i + 1]);
  }

  return result;
}
