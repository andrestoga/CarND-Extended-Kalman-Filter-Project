#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

// #define DEBUG_DATA 

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
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

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
   * Finish initializing the FusionEKF.
   * Set the process and measurement noises
   */

  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  ekf_.I_ = MatrixXd::Identity(4, 4);

  // state transition matrix
  ekf_.F_ = MatrixXd(4, 4);

  // state covariance matrix
  ekf_.P_ = MatrixXd(4, 4);

  // process covariance matrix
  ekf_.Q_ = MatrixXd(4, 4);
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */

  static int ekf_iteration = 0;

   VectorXd z = measurement_pack.raw_measurements_;

   #ifdef DEBUG_DATA

   cout << "EKF iteration: " << ekf_iteration << endl; 
   cout << "Data from: " << measurement_pack.sensor_type_ << " sensor" << endl; 
   cout << z << endl;
   cout << endl;

   #endif

  ekf_iteration++;

  if (!is_initialized_) {
    /**
     * Initialize the state ekf_.x_ with the first measurement.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    double x = 0.0;
    double y = 0.0;

    previous_timestamp_ = measurement_pack.timestamp_;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates
      //         and initialize state.
     double rho = z(0);
     double theta = z(1);
     x = sin(theta) * rho;
     y = cos(theta) * rho;
     // TODO: Check if I need to set values in the velocity components
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // Initialize state.
     x = z(0);
     y = z(1);
    }

     ekf_.x_(0) = x;
     ekf_.x_(1) = y;

     // * Create the covariance matrix.
    ekf_.P_ << 1, 0,    0,    0,
               0, 1,    0,    0,
               0, 0, 1000,    0,
               0, 0,    0, 1000;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }


  /**
   * Prediction
   */

  /**
   * Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  // Converting the timestamp from us to s
  float delta_T = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  float delta_T2 = delta_T * delta_T;
  float delta_T3 = delta_T2 * delta_T;
  float delta_T4 = delta_T3 * delta_T;

  ekf_.F_ << 1, 0, delta_T, 0,
             0, 1, 0, delta_T,
             0, 0, 1, 0,
             0, 0, 0, 1;

  double noise_ax = 9.0;
  double noise_ay = 9.0;

  ekf_.Q_ << delta_T4/4.0*noise_ax, 0.0, delta_T3/2.0*noise_ax, 0.0,
             0.0, delta_T4/4.0*noise_ay, 0.0, delta_T3/2.0*noise_ay,
             delta_T3/2.0*noise_ax, 0.0, delta_T2*noise_ax, 0.0,
             0.0, delta_T3/2.0*noise_ay, 0.0, delta_T2*noise_ay;

  ekf_.Predict();

   // #ifndef DEBUG_DATA
  /**
   * Update
   */

  /**
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    try
    {
      // Radar updates
      ekf_.R_ = R_radar_;
      ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
      ekf_.UpdateEKF(z);
    }
    catch(const std::invalid_argument&)
    {
      // Comment out this 
      cout << "Skipping this update due to division by 0 in the jacobian" << endl;
    }

  } else {
    // Laser updates
    ekf_.R_ = R_laser_;
    ekf_.H_ = H_laser_; 
    ekf_.Update(z);
  }

  // #endif

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;


  // getchar();
}
