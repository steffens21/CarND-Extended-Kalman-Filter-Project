#include <iostream>
#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  // measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0,      0.0225;

  // measurement covariance matrix - radar
  R_radar_ << 0.09, 0,      0,
              0,    0.0009, 0,
              0,    0,      0.09;

  /**
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */

  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  // Create a 4D state vector, we don't know yet the values of the x state
  Eigen::VectorXd x_in = VectorXd(4);

  // State covariance matrix P
  Eigen::MatrixXd P_in = MatrixXd(4, 4);
  P_in << 1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1000, 0,
          0, 0, 0, 1000;

  // Measurement covariance
  Eigen::MatrixXd R_in = MatrixXd(2, 2);

  // Measurement matrix
  Eigen::MatrixXd H_in = MatrixXd(2, 4);

  // The initial transition matrix F_
  Eigen::MatrixXd F_in = MatrixXd(4, 4);
  F_in << 1, 0, 1, 0,
          0, 1, 0, 1,
          0, 0, 1, 0,
          0, 0, 0, 1;

  Eigen::MatrixXd Q_in = MatrixXd(4, 4);

  ekf_.Init(x_in, P_in, F_in, H_in, R_in, Q_in);

  // Set the acceleration noise components
  noise_ax = 9;
  noise_ay = 9;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if ( !is_initialized_ ) {
    /**
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // First measurement
    std::cout << "EKF: " << std::endl;

    if ( measurement_pack.sensor_type_ == MeasurementPackage::RADAR ) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      float rho_dot = measurement_pack.raw_measurements_[2];
      float px = rho * cos(phi);
      float py = rho * sin(phi);
      // The last 2 components can also be initalized to 0, 0
      ekf_.x_ << px, py, rho_dot*cos(phi), rho_dot*sin(phi);
    } else if ( measurement_pack.sensor_type_ == MeasurementPackage::LASER ) {
      /**
      Initialize state.
      */
      ekf_.x_ << measurement_pack.raw_measurements_[0],
                 measurement_pack.raw_measurements_[1],
                 0, 0;
    }

    previous_timestamp_ = measurement_pack.timestamp_;

    // Done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  // Compute the time elapsed between the current and previous measurements
  // dt is in seconds
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  // Set the process covariance matrix Q
  ekf_.Q_ <<  dt_4/4*noise_ax,    0,               dt_3/2*noise_ax, 0,
              0,                  dt_4/4*noise_ay, 0,               dt_3/2*noise_ay,
              dt_3/2*noise_ax,    0,               dt_2*noise_ax,   0,
              0,                  dt_3/2*noise_ay, 0,               dt_2*noise_ay;

  // Modify the F matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */
  if ( measurement_pack.sensor_type_ == MeasurementPackage::RADAR ) {
    // Radar updates
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.R_ = R_laser_;
    ekf_.H_ = H_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // Print the output
  std::cout << "x_ = " << ekf_.x_ << std::endl;
  std::cout << "P_ = " << ekf_.P_ << std::endl;
}
