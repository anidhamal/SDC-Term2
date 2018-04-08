#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>       /* atan2 */

#define PI 3.14159265
using namespace std;


using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);


  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_ = false;
  previous_timestamp_ = 0.0;

  ///* State dimension
  n_x_ = 5;

  ///* Augmented state dimension
  n_aug_ = 7;

  ///* Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  MatrixXd R_laser_;// = MatrixXd( 2, 2);

  MatrixXd R_radar_;// = MatrixXd(3, 3);

  /*F_ << 1, 0, 1, 0,
  		0, 1, 0, 1,
  		0, 0, 1, 0,
  		0, 0, 0, 1;

  Q_ = MatrixXd(4, 4);*/

  P_ = MatrixXd(5, 5);
  P_ <<  1, 0, 0, 0, 0,
  		 0, 1, 0, 0, 0,
  		 0, 0, 1, 0, 0,
  		 0, 0, 0, 1, 0,
  		 0, 0, 0, 0, 1;

  weights_ = VectorXd(2*n_aug_+1);
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  std::cout << "********************************************" << std::endl << x_ << std::endl;
  std::cout << "Sensor Type : " << std::endl << measurement_pack.sensor_type_ << std::endl;
  cout << "raw_measurements_: " << measurement_pack.raw_measurements_ << endl;
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    cout << "initializing: "<< endl;
    // first measurement
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double rho = measurement_pack.raw_measurements_[0]; // range
      double phi = measurement_pack.raw_measurements_[1]; // bearing
      double rho_dot = measurement_pack.raw_measurements_[2]; // velocity of rh

      double vx = rho_dot * cos(phi);
      double vy = rho_dot * sin(phi);
      double v = sqrt(vx * vx + vy * vy);
      x_ <<  rho * cos(phi), // px
             rho * sin(phi), // py
                 v         , // v
                 0         , // yaw angle
                 0         ; // yaw rate
      cout << "RADAR: " << endl;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_ << measurement_pack.raw_measurements_[0], // px
            measurement_pack.raw_measurements_[1], // py
                            0                     , // v
                            0                     , // yaw angle
                            0                     ; // yaw rate
      cout << "LASER: " << endl;
    }
    previous_timestamp_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    cout << "x: " << x_ << endl;
    cout << "initialized: " << endl;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;
  Prediction(dt);

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    UpdateRadar(measurement_pack);
  } else {
    UpdateLidar(measurement_pack);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  **/
  AugmentedSigmaPoints(&Xsig_aug_);
  SigmaPointPrediction(&Xsig_pred_, delta_t);
  PredictMeanAndCovariance(&x_, &P_);

  cout << "x_: " << x_ << endl;
  cout << "P_: " << P_ << endl;
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

  std::cout << "UpdateLidar: "<< std::endl;
  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 2;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //transform sigma points into measurement space
  //calculate mean predicted measurement
  //calculate innovation covariance matrix S
  //MatrixXd Z_pred = MatrixXd(n_z, 2 * n_aug_ + 1);
  z_pred.fill(0);
  for (int i =  0; i< 2*n_aug_+1; i++) {
      double  px = Xsig_pred_(0,i);
      double  py = Xsig_pred_(1,i);

      VectorXd x = VectorXd(n_z);
               x << px, py;
      Zsig.col(i) = x;
      z_pred = z_pred + x * weights_(i);
  }

  MatrixXd R = MatrixXd(2, 2);
           R << std_laspx_ * std_laspx_,            0            ,
                           0           , std_laspy_ * std_laspy_ ;

  std::cout << "R: "<< R << std::endl;
  S.fill(0);
  //std::cout << "S: "<< S << std::endl;
  for (int i =0; i< 2*n_aug_+1; i++) {
      VectorXd diff = (Zsig.col(i) - z_pred);

      //std::cout << "diff: "<< diff << std::endl;
      S = S + weights_(i) * diff * diff.transpose();
  }

  S = S + R;
/*******************************************************************************
 * Student part end
 ******************************************************************************/

  //print result
  std::cout << "z_pred: " << std::endl << z_pred << std::endl;
  std::cout << "S: " << std::endl << S << std::endl;

  //write result
//  *z_out = z_pred;
//  *S_out = S;


  //create example vector for incoming radar measurement
  VectorXd z = meas_package.raw_measurements_;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //calculate cross correlation matrix
  Tc.fill(0);
  for (int i=0; i < 2 * n_aug_ +1; i++) {
      MatrixXd temp = MatrixXd(n_x_, n_z);
      VectorXd sig_diff = Xsig_pred_.col(i) - x_;
      VectorXd pred_diff = Zsig.col(i) - z_pred;
      temp = weights_(i) * sig_diff * pred_diff.transpose();
      Tc += temp;
  }
  //calculate Kalman gain K;
  MatrixXd K = MatrixXd(n_x_, n_z);
  K = Tc *  S.inverse();

  VectorXd z_diff = (z - z_pred);
  VectorXd NIS_lidar = z_diff.transpose() * S.inverse() * z_diff;
  //update state mean and covariance matrix
  x_ = x_ + K * (z - z_pred);
  P_ =  P_ - K * S * K.transpose();

/*******************************************************************************
 * Student part end
 ******************************************************************************/

  //print result
  std::cout << "Updated state x: " << std::endl << x_ << std::endl;
  std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;
  std::cout << "Normalized Innovation Squared (lidar): " << std::endl << NIS_lidar << std::endl;

  //write result
  //*x_out = x;
  //*P_out = P;
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
    //set measurement dimension, radar can measure r, phi, and r_dot
    int n_z = 3;

    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);

    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);

  /*******************************************************************************
   * Student part begin
   ******************************************************************************/

    //transform sigma points into measurement space
    //calculate mean predicted measurement
    //calculate innovation covariance matrix S
    //MatrixXd Z_pred = MatrixXd(n_z, 2 * n_aug_ + 1);
    z_pred.fill(0);
    for (int i =  0; i< 2*n_aug_+1; i++) {
        double  px = Xsig_pred_(0,i);
        double  py = Xsig_pred_(1,i);
        double  v = Xsig_pred_(2,i);
        double  psi = Xsig_pred_(3,i);

        double rho = sqrt(px * px + py *py);
        double phi = atan2(py, px);
        double rho_dot = (px * cos(psi) * v + py * sin(psi) * v)/ rho;
        VectorXd x = VectorXd(n_z);
                 x << rho, phi, rho_dot;
        Zsig.col(i) = x;
        z_pred = z_pred + x * weights_(i);
    }

    MatrixXd R = MatrixXd(3, 3);
             R << std_radr_ * std_radr_,            0            ,           0          ,
                            0          ,std_radphi_ * std_radphi_,           0          ,
                            0          ,            0            ,std_radrd_ * std_radrd_;
    //std::cout << "S: "<< S << std::endl;
    S.fill(0);
    //std::cout << "S: "<< S << std::endl;
    for (int i =0; i< 2*n_aug_+1; i++) {
        VectorXd diff = (Zsig.col(i) - z_pred);
        while(diff[1]> PI || diff[1] < -PI) {
            if (diff[1] > PI) {
                diff[1] -= 2*PI;
            } else {
                diff[1] += 2*PI;
            }
        }
        //std::cout << "diff: "<< diff << std::endl;
        S = S + weights_(i) * diff * diff.transpose();
    }

    S = S + R;
  /*******************************************************************************
   * Student part end
   ******************************************************************************/

    //print result
    std::cout << "z_pred: " << std::endl << z_pred << std::endl;
    std::cout << "S: " << std::endl << S << std::endl;

    //write result
  //  *z_out = z_pred;
  //  *S_out = S;


    //create example vector for incoming radar measurement
    VectorXd z = meas_package.raw_measurements_;

    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);

  /*******************************************************************************
   * Student part begin
   ******************************************************************************/

    //calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i=0; i < 2 * n_aug_ +1; i++) {
        MatrixXd temp = MatrixXd(n_x_, n_z);
        VectorXd sig_diff = Xsig_pred_.col(i) - x_;
        //angle normalization
        while (sig_diff(3)> PI) sig_diff(3)-= 2.* PI;
        while (sig_diff(3)< -PI) sig_diff(3)+= 2.* PI;

        VectorXd pred_diff = Zsig.col(i) - z_pred;
        while (pred_diff(1)> PI) pred_diff(1)-= 2.*PI;
        while (pred_diff(1)<-PI) pred_diff(1)+= 2.*PI;

        temp = weights_(i) * sig_diff * pred_diff.transpose();
        Tc += temp;
    }
    //calculate Kalman gain K;
    MatrixXd K = MatrixXd(n_x_, n_z);
    K = Tc *  S.inverse();

    VectorXd z_diff = z - z_pred;
    //angle normalization
    while (z_diff(1)> PI) z_diff(1) -= 2.*PI;
    while (z_diff(1)<-PI) z_diff(1) += 2.*PI;

    //calculate NIS
    VectorXd NIS_radar = z_diff.transpose() * S.inverse() * z_diff;
    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ =  P_ - K * S * K.transpose();

  /*******************************************************************************
   * Student part end
   ******************************************************************************/

    //print result
    std::cout << "Updated state x: " << std::endl << x_ << std::endl;
    std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;
    std::cout << "Normalized Innovation Squared (radar): " << std::endl << NIS_radar << std::endl;
}

void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug = MatrixXd::Zero(n_aug_, n_aug_);

  //create sigma point matrix
  //MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //create augmented mean state
  x_aug.head(n_x_) = x_;
  x_aug[n_x_] = 0;
  x_aug[n_x_+1] = 0;
  //std::cout << "x_aug = " << std::endl << x_aug << std::endl;
  //create augmented covariance matrix
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  MatrixXd Q = MatrixXd(2,2);
           Q << std_a_ * std_a_,          0,
                      0      , std_yawdd_ * std_yawdd_;
  P_aug.bottomRightCorner(2, 2) = Q;
  //P_aug(n_x+1, n_x+1) = std_yawdd * std_yawdd;
  //std::cout << "P_aug = " << std::endl << P_aug << std::endl;
  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug_.col(0) = x_aug;
  for (int i=0; i<n_aug_; i++) {
      VectorXd l = (sqrt(lambda_ + n_aug_)) * A.col(i);
      Xsig_aug_.col(i+1) = x_aug + l;
      Xsig_aug_.col(i+1+n_aug_) = x_aug - l;
  }
/*******************************************************************************
 * Student part end
 ******************************************************************************/

  //print result
  std::cout << "Xsig_aug = " << std::endl << Xsig_aug_ << std::endl;

  //write result
  //*Xsig_out = Xsig_aug;
}

void UKF::SigmaPointPrediction(MatrixXd* Xsig_out, double dt) {
  //create matrix with predicted sigma points as columns
  //MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

  double delta_t = dt; //time diff in sec
/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //predict sigma points
  //avoid division by zero
  //write predicted sigma points into right column
  for (int i=0; i < (2 * n_aug_ + 1); i++) {
      VectorXd x = VectorXd(5);
      double    p_x  = Xsig_aug_(0, i);
      double    p_y  = Xsig_aug_(1, i);
      double     v   = Xsig_aug_(2, i);
      double    psi  = Xsig_aug_(3, i);
      double psi_dot = Xsig_aug_(4, i);
      double    v_a  = Xsig_aug_(5, i);
      double vpsi_d_2= Xsig_aug_(6, i);

      double r0;
      double r1;
      if (psi_dot != 0) {
          r0 = (v/ psi_dot )* (sin(psi + psi_dot *delta_t) - sin(psi)) + 0.5 * delta_t * delta_t * cos(psi) * v_a;
          r1 = (v/ psi_dot )* (cos(psi) - cos(psi + psi_dot *delta_t)) + 0.5 * delta_t * delta_t * sin(psi) * v_a;
      } else {
          r0 = v * cos(psi) * delta_t + 0.5 * delta_t * delta_t * cos(psi) * v_a;
          r1 = v * sin(psi) * delta_t + 0.5 * delta_t * delta_t * sin(psi) * v_a;
      }
      double r2 = delta_t * v_a;
      double r3 = psi_dot * delta_t + 0.5 * delta_t * delta_t * vpsi_d_2;
      double r4 = delta_t * vpsi_d_2;
      x << p_x + r0, p_y + r1, v + r2, psi + r3, psi_dot + r4;
      Xsig_pred_.col(i) = x;
  }

/*******************************************************************************
 * Student part end
 ******************************************************************************/

  //print result
  std::cout << "Xsig_pred = " << std::endl << Xsig_pred_ << std::endl;

  //write result
  //*Xsig_out = Xsig_pred;
}

void UKF::PredictMeanAndCovariance(VectorXd* x_out, MatrixXd* P_out) {
 std::cout << "PredictMeanAndCovariance = " << std::endl;
  //create vector for predicted state
  //VectorXd x = VectorXd(n_x_);

  //create covariance matrix for prediction
  //MatrixXd P = MatrixXd(n_x_, n_x_);


/*******************************************************************************
 * Student part begin
 ******************************************************************************/
  // set weights
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }
  std::cout << "x = "<< x_ << std::endl;
  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    std::cout << "x_diff = " << x_diff(3) << std::endl;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-= 2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+= 2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }


/*******************************************************************************
 * Student part end
 ******************************************************************************/

  //print result
  std::cout << "Predicted state" << std::endl;
  std::cout << x_ << std::endl;
  std::cout << "Predicted covariance matrix" << std::endl;
  std::cout << P_ << std::endl;

  //write result
  //*x_out = x;
  //*P_out = P;
}