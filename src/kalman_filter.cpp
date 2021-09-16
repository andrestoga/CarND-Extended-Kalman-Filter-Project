#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/*
 * Please note that the Eigen library does not initialize
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
   * predict the state
   */
   x_ = F_ * x_;
   P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * update the state by using Kalman Filter equations
   */
   MatrixXd y = z - H_ * x_; // 2x1 - 2x4 * 4x1 = 2x1
   MatrixXd HT = H_.transpose();
   MatrixXd S = H_ * P_ * HT + R_; // 2x4 * 4x4 * 4x2 + 2x2 = 2x2
   MatrixXd K = P_ * HT * S.inverse(); // 4x4 * 4x2 * 2x2 = 4x2

   x_ = x_ + K * y; // 4x1 + 4x2 * 2x1 = 4x1
   P_ = (I_ - K * H_) * P_; // (4x4 - 4x2 * 2x4) * 4x4 = 4x4
}

// void KalmanFilter::Update(const MatrixXd& H)
// {
//    MatrixXd y = z - H * x_;
//    MatrixXd S = H * P_ * H.transpose() + R_;
//    MatrixXd K = P_ * H.transpose() * S.transpose();
//
//    x_ = x_ + K * y;
//    P_ = (I_ - K * H) * P_;
// }

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * update the state by using Extended Kalman Filter equations
   */
   MatrixXd y = z - h(x_); // 3x1 - 3x1 = 3x1

   // Check the angle that it is in the interval -pi to pi 
   double angle = y(1);

   if (angle > M_PI)
   {
     angle -= 2.0 * M_PI;
   }
   else if (angle < -M_PI)
   {
     angle += 2.0 * M_PI;
   }

   y(1) = angle;

   MatrixXd HT = H_.transpose();
   MatrixXd S = H_ * P_ * HT + R_; // 3x4 * 4x4 * 4x3 + 3x3 = 3x3
   MatrixXd K = P_ * HT * S.inverse(); // 4x4 * 4x3 * 3x3 = 4x3

   x_ = x_ + K * y; // 4x1 + 4x3 * 3x1 = 4x1
   P_ = (I_ - K * H_) * P_; // (4x4 - 4x3 * 3x4) * 4x4 = 4x4
}

VectorXd KalmanFilter::h(const VectorXd x)
{
  VectorXd res(3);

  double px = x(0);
  double py = x(1);
  double vx = x(2);
  double vy = x(3);

  double px2 = px*px;
  double py2 = py*py;

  double c = sqrt(px2 + py2);

  res(0) = c; 
  res(1) = atan2(py, px); 
  res(2) = (px * vx + py * vy) / c;

  return res; 
}