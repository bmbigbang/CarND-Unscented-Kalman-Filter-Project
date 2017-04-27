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
  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.35;

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

  // set the linear laser measurement matrix
  H_laser_ << 1, 0, 0, 0, 0,
              0, 1, 0, 0, 0;

  // set the linear radar measurement matrix
  H_radar_ << 0, 0, 1, 0, 0,
              0, 0, 0, 1, 0,
              0, 0, 0, 0, 1;

  // set state dimension
  n_x_ = 5;

  // set augmented dimension
  n_aug_ = 7;

  // set measurement dimension, radar can measure r, phi, and r_dot
  n_z_ = 3;

  // define spreading parameter
  lambda_ = 3 - n_aug_;

  //create augmented mean vector
  x_aug_ = VectorXd(7);

  //create augmented state covariance
  P_aug_ = MatrixXd(7, 7);

  //create sigma point matrix
  Xsig_ = MatrixXd(n_x_, 2 * n_x_ + 1);

  //create sigma point matrix
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
  Zsig_ = MatrixXd(n_z_, 2 * n_aug_ + 1);

  // mean predicted measurement
  z_pred_ = VectorXd(n_z_);
  z_pred_.fill(0.0);

  ///* the current NIS for radar and 0.5 percent ki threshold
  double NIS_radar_;
  NIS_r_tresh_ = 7.815;

  ///* the current NIS for laser and 0.5 percent ki threshold
  double NIS_laser_;
  NIS_l_tresh_ = 5.991;

  // process noise covariance matrix
  MatrixXd Q_ = MatrixXd(5, 5);

  // initialize state transition matrix
  MatrixXd F_ = MatrixXd(5, 5);
  F_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;

  // identity matrix for calculations
  MatrixXd I_ = MatrixXd::Identity(5, 5);

  // initialize place holders for kalman filter calculation
  S_ = MatrixXd(5, 5);
  S_.fill(0.0);

  K_ = MatrixXd(5, 5);
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
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      //set the state with the initial location and zero velocity and yaw

      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0.0, 0.0, 0.0;
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
  }


  /*****************************************************************************
 *  Update
 ****************************************************************************/

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates

    float px = x_[0];
    float py = x_[1];
    float v = x_[2];
    float psi = x_[3];
    float psidot = x_[4];

    //check division by zero
    if (px <= 0.00001) {
      px = 0.00001;
    }
    if (py <= 0.00001) {
      py = 0.00001;
    }

    float rho = sqrt(px * px + py * py);
    float phi = atan(py / px);
    float rho_dot = (px * vx + py * vy) / rho;

    ekf_.m_ << rho, phi, rho_dot;

    UpdateRadar(meas_package);
  } else {
    // Laser updates
    UpdateLidar(meas_package);
  }

  // print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
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

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig_(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig_(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig_(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  for (int i=0; i < 2*n_aug_+1; i++) {
    z_pred_ = z_pred_ + weights_(i) * Zsig_.col(i);
  }

  //measurement covariance matrix S
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig_.col(i) - z_pred_;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S_ = S_ + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  S_ = S_ + R_radar_;
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
}
