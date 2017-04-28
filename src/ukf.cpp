#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initialize values for the object covariance matrix
  P_ = MatrixXd(5, 5);
  P_ << 1000, 0, 0, 0, 0,
        0, 10, 0, 0, 0,
        0, 0, 10, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1000;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2.2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.22;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  // initializing matrices noise matrices and The State-Covariance translation Matrices H
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 5);
  H_radar_ = MatrixXd(3, 5);

  // measurement covariance matrix - laser
  R_laser_ << pow(std_laspx_, 2), 0,
              0, pow(std_laspy_, 2);

  // measurement covariance matrix - radar
  R_radar_ << pow(std_radr_, 2), 0, 0,
              0, pow(std_radphi_, 2), 0,
              0, 0, pow(std_radrd_, 2);

  // set state dimension
  n_x_ = 5;

  // set augmented dimension
  n_aug_ = 7;

  // set measurement dimension, radar can measure r, phi, and r_dot
  n_z_r_ = 3;
  n_z_l_ = 2;

  // define spreading parameter
  lambda_ = 3 - n_aug_;

  // create augmented mean vector
  x_aug_ = VectorXd(7);

  // create x predicted sigma points
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // create augmented state covariance
  P_aug_ = MatrixXd(7, 7);

  // create sigma point matrix
  Xsig_ = MatrixXd(n_x_, 2 * n_x_ + 1);

  // create sigma point matrix
  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // set vector for weights
  weights_ = VectorXd(2*n_aug_+1);
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  // create matrix for sigma points in measurement space
  Zsig_r_ = MatrixXd(n_z_r_, 2 * n_aug_ + 1);
  Zsig_r_.fill(0.0);
  Zsig_l_ = MatrixXd(n_z_l_, 2 * n_aug_ + 1);
  Zsig_l_.fill(0.0);

  // mean predicted measurement
  z_pred_r_ = VectorXd(n_z_r_);
  z_pred_l_ = VectorXd(n_z_l_);

  ///* the current NIS for radar and 5 percent ki threshold
  NIS_r_tresh_ = 7.815;

  ///* the current NIS for laser and 5 percent ki threshold
  NIS_l_tresh_ = 5.991;

  // identity matrix for calculations
  MatrixXd I_ = MatrixXd::Identity(5, 5);

  // initialize place holders for kalman filter calculation
  S_r_ = MatrixXd(n_z_r_, n_z_r_);
  S_r_.fill(0.0);

  S_l_ = MatrixXd(n_z_l_, n_z_l_);
  S_l_.fill(0.0);

  // cross correlation tensors
  Tc_r_ = MatrixXd(n_x_, n_z_r_);
  Tc_l_ = MatrixXd(n_x_, n_z_l_);

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {

    // first measurement
    cout << "UKF: " << endl;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      //set the state with the initial location and velocity and yaw

      double ro = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double ro_dot = meas_package.raw_measurements_[2];
      x_ << ro * cos(phi), ro * sin(phi), ro_dot, phi, 0.0;
      z_pred_r_ << ro, phi, ro_dot;

    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      //set the state with the initial location and zero velocity and yaw

      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0.0, 0.0, 0.0;
      z_pred_l_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];
    }

    previous_timestamp_ = meas_package.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }
  //dt - expressed in seconds
  double dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = meas_package.timestamp_;

  //calculate square root of P
  MatrixXd A = P_.llt().matrixL();

  //set first column of sigma point matrix
  Xsig_.col(0)  = x_;

  //set remaining sigma points
  for (int i = 0; i < n_x_; i++)
  {
    Xsig_.col(i+1)     = x_ + sqrt(lambda_+n_x_) * A.col(i);
    Xsig_.col(i+1+n_x_) = x_ - sqrt(lambda_+n_x_) * A.col(i);
  }

  //create augmented mean state
  x_aug_.head(n_x_) = x_;
  x_aug_(5) = 0;
  x_aug_(6) = 0;

  //create augmented covariance matrix
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(n_x_, n_x_) = P_;
  P_aug_(n_x_, n_x_) = pow(std_a_, 2);
  P_aug_(n_x_ + 1, n_x_ + 1) = pow(std_yawdd_, 2);

  //create square root matrix
  MatrixXd AA = P_aug_.llt().matrixL();

  //create augmented sigma points
  Xsig_aug_.col(0)  = x_aug_;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug_.col(i+1)       = x_aug_ + sqrt(lambda_+n_aug_) * AA.col(i);
    Xsig_aug_.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_+n_aug_) * AA.col(i);
  }

  // avoid predicting twice if the measurements coincide
  if ( dt > 0.001 )
  {
    Prediction(dt);
    //transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

      // extract values for better readibility
      double p_x = Xsig_pred_(0, i);
      double p_y = Xsig_pred_(1, i);
      double v = Xsig_pred_(2, i);
      double yaw = Xsig_pred_(3, i);

      double v1 = cos(yaw) * v;
      double v2 = sin(yaw) * v;

      //check division by zero
      if (p_x <= 0.00001) {
        p_x = 0.00001;
      }
      if (p_y <= 0.00001) {
        p_y = 0.00001;
      }

      // measurement model
      Zsig_r_(0,i) = sqrt(p_x*p_x + p_y*p_y);                        // p_x
      Zsig_r_(1,i) = atan2(p_y,p_x);                                 //p_y
      Zsig_r_(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
      // laser
      Zsig_l_(0, i) = p_x;                      //px
      Zsig_l_(1, i) = p_y;                       //py
    }

    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {

      //mean predicted measurement
      for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        z_pred_l_ = z_pred_l_ + weights_(i) * Zsig_l_.col(i);
      }

      //measurement covariance matrix S
      for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        //residual
        VectorXd z_diff = Zsig_l_.col(i) - z_pred_l_;

        S_l_ = S_l_ + weights_(i) * z_diff * z_diff.transpose();
      }

      S_l_ = S_l_ + R_laser_;

    } else {

      //mean predicted measurement
      for (int i=0; i < 2*n_aug_+1; i++) {
        z_pred_r_ = z_pred_r_ + weights_(i) * Zsig_r_.col(i);
      }

      //measurement covariance matrix S
      for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        //residual
        VectorXd z_diff = Zsig_r_.col(i) - z_pred_r_;

        //angle normalization
        while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
        while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

        S_r_ = S_r_ + weights_(i) * z_diff * z_diff.transpose();
      }

      S_r_ = S_r_ + R_radar_;
    }
  }


  /*****************************************************************************
 *  Update
 ****************************************************************************/

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    UpdateRadar(meas_package);
  } else {
    // Laser updates
    UpdateLidar(meas_package);
  }

  // print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
  cout << "NIS L,R = " << NIS_laser_ << "," << NIS_radar_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  for (int i = 0; i< (2*n_aug_) +1; i++) {
    //extract values for better readability
    double p_x = Xsig_aug_(0, i);
    double p_y = Xsig_aug_(1, i);
    double v = Xsig_aug_(2, i);
    double yaw = Xsig_aug_(3, i);
    double yawd = Xsig_aug_(4, i);
    double nu_a = Xsig_aug_(5, i);
    double nu_yawdd = Xsig_aug_(6, i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_p = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
      py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
    } else {
      px_p = p_x + v * delta_t * cos(yaw);
      py_p = p_y + v * delta_t * sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_p = v_p + nu_a * delta_t;

    yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p = yawd_p + nu_yawdd * delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;
  }

  //predict state mean
  for (int i = 0; i < (2 * n_aug_) + 1; i++) {
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  //predict state covariance matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) {x_diff(3)-=2.*M_PI;}
    while (x_diff(3)<-M_PI) {x_diff(3)+=2.*M_PI;}

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  // calculate cross correlation matrix
  Tc_l_.fill(0.0);

  // gather measurement values
  double p_x = meas_package.raw_measurements_[0];
  double p_y = meas_package.raw_measurements_[1];
  z_ = VectorXd(2);
  z_ << p_x, p_y;

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig_l_.col(i) - z_pred_l_;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc_l_ = Tc_l_ + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  K_ = Tc_l_ * S_l_.inverse();

  //residual
  VectorXd z_diff = z_ - z_pred_l_;

  // calculate NIS
  NIS_laser_ = z_diff.transpose() * S_l_.inverse() * z_diff;

  //update state mean and covariance matrix
  x_ = x_ + K_ * z_diff;
  P_ = P_ - K_*S_l_*K_.transpose();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  // calculate cross correlation matrix
  Tc_r_.fill(0.0);
  cout << meas_package.raw_measurements_ << endl;
  // gather measurement values
  double ro = meas_package.raw_measurements_[0];
  double phi = meas_package.raw_measurements_[1];
  double ro_dot = meas_package.raw_measurements_[2];
  z_ = VectorXd(3);
  z_ << ro, ro_dot, phi;

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig_r_.col(i) - z_pred_r_;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc_r_ = Tc_r_ + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  K_ = Tc_r_ * S_r_.inverse();

  //residual
  VectorXd z_diff = z_ - z_pred_r_;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  // calculate NIS
  NIS_radar_ = z_diff.transpose() * S_r_.inverse() * z_diff;

  //update state mean and covariance matrix
  x_ = x_ + K_ * z_diff;
  P_ = P_ - K_*S_r_*K_.transpose();
}
