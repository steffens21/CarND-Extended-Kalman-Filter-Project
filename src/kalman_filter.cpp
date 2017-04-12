#include "kalman_filter.h"
#include <math.h>
#include "tools.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

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
  TODO:
    * predict the state
  */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_; // TODO: use R_laser here
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */

  float rho = sqrt(x_[0]*x_[0] + x_[1]*x_[1]);
  float phi = atan2(x_[1], x_[0]);
  float rhodot = (x_[0]*x_[2] + x_[1]*x_[3]) / rho;
  //cout << "rho: " << rho << "\tphi: " << phi << "\trhodot: " << rhodot << endl;
  VectorXd z_pred(3);
  z_pred << rho , phi, rhodot;
  VectorXd y = z - z_pred;
  //cout << "y: " << y << endl;
  Tools tool;
  MatrixXd Hj = tool.CalculateJacobian(x_);
  //cout << "Hj: " << Hj << endl;
  MatrixXd Hjt = Hj.transpose();
  MatrixXd S = Hj * P_ * Hjt + R_;
  //cout << "S: " << S << endl;
  Eigen::FullPivLU<MatrixXd> lu(S);
  if (!lu.isInvertible())
  {
    return void();
  }
  MatrixXd Si = S.inverse();
  //cout << "Si: " << Si << endl;
  MatrixXd PHjt = P_ * Hjt;
  MatrixXd K = PHjt * Si;
  //cout << "K: " << K << endl;

  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * Hj) * P_;
}
