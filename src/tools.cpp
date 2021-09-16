#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * Calculate the RMSE here.
   */
  VectorXd res(4);

  res << 0, 0, 0, 0;
  
  if (estimations.size() != ground_truth.size())
  {
   return res; 
  }

  for (int i = 0; i < estimations.size(); ++i)
  {
    VectorXd tmp = estimations[i].array() - ground_truth[i].array();
    tmp = tmp.array() * tmp.array();

    res += tmp;
  }

  res = res / estimations.size();
  res = res.array().sqrt();

  return res;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * Calculate a Jacobian here.
   */
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);

  MatrixXd jaco = MatrixXd(3,4);

  jaco << 0, 0, 0, 0,
          0, 0, 0, 0,
          0, 0, 0, 0;
  
  double px2 = px * px;
  double py2 = py * py;
  double px2_py2 = px2 + py2;
  double c1 = sqrt(px2_py2);
  double c2 = pow(px2_py2, 3/2);

  // TODO: I don't know what value to assign to the jacobian
  // when there is division by 0
  if (px2_py2 > 0.0001)
  {
    jaco << px/c1, py/c1, 0, 0,
            -py/px2_py2, px/px2_py2, 0, 0,
            py*(vx*py - vy*px)/c2, px*(vy*px - vx*py)/c2, px/c1, py/c2;
  }
  else
  {
   throw std::invalid_argument("Jacobian division by 0");
  }

  return jaco;
}
