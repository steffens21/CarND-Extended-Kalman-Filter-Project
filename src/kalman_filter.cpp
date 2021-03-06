#include <math.h>
#include <iostream>
#include "kalman_filter.h"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

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
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  // Check if S is invertable, if not, don't update
  Eigen::FullPivLU<MatrixXd> lu(S);
  if (!lu.isInvertible()) {
    return void();
  }
  MatrixXd Si = S.inverse();
  MatrixXd K = P_ * Ht * Si;

  // new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
    * update the state by using Extended Kalman Filter equations
  */

  // Convert to polar coordinates
  float rho = sqrt(x_[0]*x_[0] + x_[1]*x_[1]);
  float phi = atan2(x_[1], x_[0]);
  float rhodot = (x_[0]*x_[2] + x_[1]*x_[3]) / rho;
  VectorXd z_pred(3);
  z_pred << rho , phi, rhodot;

  VectorXd y = z - z_pred;

  Tools tool;
  float c1 = x_(0)*x_(0)+x_(1)*x_(1);
  // Check division by zero would occur in Jacobian calculation.
  // If so, don't update x and P
  if ( fabs(c1) < 0.001 ) {
    return void();
  }
  MatrixXd Hj = tool.CalculateJacobian(x_);
  MatrixXd Hjt = Hj.transpose();
  MatrixXd S = Hj * P_ * Hjt + R_;

  // Check if S is invertable, if not, don't update
  Eigen::FullPivLU<MatrixXd> lu(S);
  if ( !lu.isInvertible() ) {
    return void();
  }
  MatrixXd Si = S.inverse();
  MatrixXd K = P_ * Hjt * Si;

  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * Hj) * P_;
}
