#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

class MPC {
 public:
  MPC();

  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Returns (delta0, a0, mpc_x, mpc_y).
  tuple<double, double, vector<double>, vector<double>> Solve(Eigen::VectorXd coeffs, double v0, double ref_v);
};

#endif /* MPC_H */
