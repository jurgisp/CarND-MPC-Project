#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

const size_t N = 20;
const double dt = 0.1;

const size_t x_start = 0;
const size_t y_start = N;
const size_t psi_start = 2 * N;
const size_t v_start = 3 * N;
const size_t cte_start = 4 * N;
const size_t epsi_start = 5 * N;
const size_t delta_start = 6 * N;
const size_t a_start = 6 * N + (N - 1);
const size_t n_vars = 6 * N + 2 * (N - 1);
const size_t n_constraints = 6 * N;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

AD<double> polyeval(Eigen::VectorXd coeffs, AD<double> x) {
  AD<double> result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * CppAD::pow(x, i);
  }
  return result;
}

AD<double> polyderiveval(Eigen::VectorXd coeffs, AD<double> x) {
  AD<double> result = 0.0;
  for (int i = 1; i < coeffs.size(); i++) {
    result += i * coeffs[i] * CppAD::pow(x, i-1);
  }
  return result;
}

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  double v0;
  double ref_v;

  FG_eval(Eigen::VectorXd coeffs, double v0, double ref_v) {
    this->coeffs = coeffs;
    this->v0 = v0;
    this->ref_v = ref_v;
  }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {

    // Cost f[0]

    fg[0] = 0;

    // The part of the cost based on the reference state.
    for (size_t t = 0; t < N; t++) {
      fg[0] += CppAD::pow(vars[cte_start + t], 2);
      fg[0] += CppAD::pow(vars[epsi_start + t], 2);
      fg[0] += CppAD::pow(vars[v_start + t] - ref_v, 2);
    }

    // Minimize the use of actuators.
    for (size_t t = 0; t < N - 1; t++) {
      fg[0] += CppAD::pow(vars[delta_start + t], 2);
      fg[0] += CppAD::pow(vars[a_start + t], 2);
    }

    // Minimize the value gap between sequential actuations.
    for (size_t t = 0; t < N - 2; t++) {
      fg[0] += CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2);
      fg[0] += CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
    }

    // Initial conditions

    fg[1 + x_start] = vars[x_start];
    fg[1 + y_start] = vars[y_start];
    fg[1 + psi_start] = vars[psi_start];
    fg[1 + v_start] = vars[v_start] - v0;
    fg[1 + cte_start] = vars[cte_start] - coeffs[0]; // coeffs[0] == value at x=0
    fg[1 + epsi_start] = vars[epsi_start] - CppAD::atan(coeffs[1]); // coeffs[1] == derivative at x=0

    // Constraints between steps

    for (size_t t = 1; t < N; t++) {
      // The state at time t+1 .
      AD<double> x1 = vars[x_start + t];
      AD<double> y1 = vars[y_start + t];
      AD<double> psi1 = vars[psi_start + t];
      AD<double> v1 = vars[v_start + t];
      AD<double> cte1 = vars[cte_start + t];
      AD<double> epsi1 = vars[epsi_start + t];

      // The state at time t.
      AD<double> x0 = vars[x_start + t - 1];
      AD<double> y0 = vars[y_start + t - 1];
      AD<double> psi0 = vars[psi_start + t - 1];
      AD<double> v0 = vars[v_start + t - 1];
      AD<double> cte0 = vars[cte_start + t - 1];
      AD<double> epsi0 = vars[epsi_start + t - 1];

      // Only consider the actuation at time t.
      AD<double> delta0 = vars[delta_start + t - 1];
      AD<double> a0 = vars[a_start + t - 1];

      AD<double> f0 = polyeval(coeffs, x0);
      AD<double> df0 = polyderiveval(coeffs, x0);
      AD<double> psides0 = CppAD::atan(df0);

      // Equations for the model:
      // x_[t] = x[t-1] + v[t-1] * cos(psi[t-1]) * dt
      // y_[t] = y[t-1] + v[t-1] * sin(psi[t-1]) * dt
      // psi_[t] = psi[t-1] + v[t-1] / Lf * delta[t-1] * dt
      // v_[t] = v[t-1] + a[t-1] * dt
      // cte[t] = f(x[t-1]) - y[t-1] + v[t-1] * sin(epsi[t-1]) * dt
      // epsi[t] = psi[t] - psides[t-1] + v[t-1] * delta[t-1] / Lf * dt
      fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[1 + psi_start + t] = psi1 - (psi0 + v0 * delta0 / Lf * dt);
      fg[1 + v_start + t] = v1 - (v0 + a0 * dt);
      fg[1 + cte_start + t] =
              cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
      fg[1 + epsi_start + t] =
              epsi1 - ((psi0 - psides0) + v0 * delta0 / Lf * dt);
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

tuple<double, double, vector<double>, vector<double>> MPC::Solve(Eigen::VectorXd coeffs, double v0, double ref_v) {
  bool ok = true;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  chrono::system_clock::time_point timer = chrono::system_clock::now();

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  // Set unlimited bounds by default
  Dvector vars(n_vars);
  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  for (size_t i = 0; i < n_vars; i++) {
    vars[i] = 0;
    vars_lowerbound[i] = -1e9;
    vars_upperbound[i] = 1e9;
  }

  // Initial state
  // Since polynomial is fit in car-relative coordinates, the initial position is by definition (0,0,0)
  vars[x_start] = 0;
  vars[y_start] = 0;
  vars[psi_start] = 0;
  vars[v_start] = v0;
  vars[cte_start] = coeffs[0]; // Value at x=0
  vars[epsi_start] = CppAD::atan(coeffs[1]); // Derivative at x=0

  // Set bounds only for actuator variables

  for (size_t t = 0; t < N-1; t++) {
    vars_lowerbound[delta_start + t] = -25 / 180.0 * M_PI;
    vars_upperbound[delta_start + t] = 25 / 180.0 * M_PI;

    vars_lowerbound[a_start + t] = -5;
    vars_upperbound[a_start + t] = 5;
  }


  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (size_t i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs, v0, ref_v);

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
  options += "Sparse  true        forward\n";
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

  if (ok) {
    auto cost = solution.obj_value;
    long elapsed_ms = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now() - timer).count();
    cout << "MPC solved in " << elapsed_ms << "ms: ";
    cout << "delta=" << solution.x[delta_start] << ", a=" << solution.x[a_start];
    cout << " cost=" << cost << endl;
  }
  else {
    // Assuming MPC succeeds, throw exception otherwise
    cout << "MPC failed - panic" << endl;
    abort();
  }

  vector<double> mpc_x;
  vector<double> mpc_y;

  for (size_t i = 0; i < N-1; i++) {
    mpc_x.push_back(solution.x[x_start+i]);
    mpc_y.push_back(solution.x[y_start+i]);
//    cout << "(x, y) = (" << solution.x[x_start+i] << ", " << solution.x[y_start+i] << ") ";
//    cout << "(v, psi) = (" << solution.x[v_start+i] << ", " << solution.x[psi_start+i] << ") ";
//    cout << "(delt, a) = (" << solution.x[delta_start+i] << ", " << solution.x[a_start+i] << ") " << endl;
  }

  return {solution.x[delta_start], solution.x[a_start], mpc_x, mpc_y};
}
