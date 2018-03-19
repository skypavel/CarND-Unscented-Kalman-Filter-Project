#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

//Special function for normalization
double angle_normalization(double angle)
{
    if (angle > M_PI){
        int kr = angle/M_PI;
        angle = angle - kr*M_PI;
        return angle;
    }

    if (angle < - M_PI){
        int kr = angle/M_PI;
        angle = angle - kr*M_PI;
        return angle;
    }
}

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

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  // init covariance matrix
  P_ = MatrixXd(5,5);
  P_ <<   1,0,0,0,0,
          0,1,0,0,0,
          0,0,1,0,0,
          0,0,0,1,0,
          0,0,0,0,1;

  // state dimension
  n_x_ = 5;

  //Augmented state dimension
  n_aug_ = 7;

  // Sigma point spreading parameter
  //lambda_ = 3 - n_aug_;
    lambda_ = 3 - n_x_;

  // Predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Create vector for weights
  weights_ = VectorXd(2 * n_aug_ + 1);

  // NIS for radar
  NIS_radar_ = 0.0;

  // NIS for lidar
  NIS_lidar_ = 0.0;

  //
  is_initialized_ = false;

  //
  time_us_ = 0.0;

  //Initialize measurement noise covariance matrix


    //add measurement noise covariance matrix for radar


}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  if (!is_initialized_) {
    //Initialize covariance matrix
      P_ << 1, 0, 0, 0, 0,
            0, 1, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;

    x_.fill(0.0);

    if (meas_package.sensor_type_ == MeasurementPackage::LASER){

      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
      double rho = meas_package.raw_measurements_(0); //range
      double phi = meas_package.raw_measurements_(1); //bearing
      double rho_dot = meas_package.raw_measurements_(2); //velocity of rho

      x_(0) = rho*cos(phi);
      x_(1) = rho*sin(phi);

    }

    // save current timestamp
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  //compute the time elapsed between the current and previous measurement
  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0; //dt - expressed in seconds
  time_us_ = meas_package.timestamp_;

  //Prediction
  Prediction(dt);

  // std::cout << "x_prediction:" << std::endl << x_ << std::endl;

  //Update
  if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    UpdateLidar(meas_package);
    std::cout << "***** Update Lidar *****"<<std::endl;
    std::cout << x_ <<std::endl;
   }
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
    UpdateRadar(meas_package);
      std::cout << "***** Update Radar *****"<<std::endl;
      std::cout << x_ <<std::endl;
  }

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  // 1. Generate Augmented Sigma Points

  //Create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //Create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

    std::cout<< "**** Matrix P in prediction ****=  "<<std::endl;
    std::cout<< P_<<std::endl;

  //Create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // Create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //Create augmented sigma points
  Xsig_aug.col(0) = x_aug;

  //Set remaining sigma points
  for (int i = 0; i < n_aug_; i++){
    Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  std::cout<< "***** yaw prediction *********"<<std::endl;

  //Predict Sigma Points
  for (int i = 0; i < 2*n_aug_ + 1; i++ ){
    //extract values for better readability
    double p_x = Xsig_aug(0, i);
    double p_y = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);


    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_p = p_x + v/yawd * (sin(yaw + yawd*delta_t) - sin(yaw));
      py_p = p_y + v/yawd * (cos(yaw) - cos(yaw + yawd*delta_t));
    }
    else {
      px_p = p_x + v*delta_t*cos(yaw);
      py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    std::cout<< i << " yaw =  "<<yaw<<" yaw_p = "<<yaw_p<<std::endl;

    //write predicted sigma points
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  //Convert Predicted Sigma Points to Mean/Covariance

  //Set weigth
  double weigth_0 = lambda_/(lambda_ + n_aug_);
  weights_(0) = weigth_0;

  for (int i = 1; i < 2*n_aug_ + 1; i++){
    double weight = 0.5/(n_aug_ + lambda_);
    weights_(i) = weight;
  }

  //Predicted state mean
  //x_.fill(0.0);
  //for (int i = 0; i < 2*n_aug_ + 1; i++){
  //  x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  //}

 std::cout << "******prediction******" << std::endl;
   x_ = Xsig_pred_ * weights_;
  while (x_(3) > M_PI) x_(3) -= 2.*M_PI;
  while (x_(3) < -M_PI) x_(3) += 2.*M_PI;
 std::cout << " x_ prediction: " << std::endl << x_ << std::endl;

  //predicted state covariance matrix
  P_.fill(0.0);
 // std::cout << "M_PI:" << M_PI << std::endl;
  for (int i = 0; i < 2*n_aug_ + 1; i++){
    //state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
   // std::cout << "x_diff:" << std::endl << x_diff << std::endl;
    while (x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;
     //x_diff(3) = angle_normalization(x_diff(3));

    P_ = P_ + weights_(i) * x_diff *x_diff.transpose();


  }
    std::cout<< "**** Matrix P after prediction ****=  "<<std::endl;
    std::cout<< P_<<std::endl;

}


/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {


  //extract measurement
  VectorXd z = meas_package.raw_measurements_;

  // set measurement dimension for lidar
  int n_z = 2;

  //create  matrix for sigma point in measurements space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //Transform sigma points into measurement space
  for (int i = 0; i < 2*n_aug_+1; i++) {
    //extract values
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);

    // measurement model
    Zsig(0, i) = p_x;
    Zsig(1, i) = p_y;
  }

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i = 0; i < 2*n_aug_ + 1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);

  for (int i = 0; i < 2*n_aug_ + 1; i++){

    VectorXd z_diff = Zsig.col(i) - z_pred;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  MatrixXd R_lidar = MatrixXd(2, 2);
  R_lidar << std_laspx_*std_laspx_, 0,
          0, std_laspy_*std_laspy_;

   S = S + R_lidar;
  std::cout<< "**** Matrix S lidar update ****=  "<<std::endl;
  std::cout<< S<<std::endl;
  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //state differnce
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    std::cout<< "**** x_diff in lidar update ****=  "<<std::endl;
    std::cout<<i<<"  "<< x_diff(3) <<" Xpred= "<< Xsig_pred_.col(i)(3)<<"  x_ = "<<x_(3)<< std::endl;
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  std::cout<< "**** Matrix Tc lidar update ****=  "<<std::endl;
  std::cout<< Tc <<std::endl;
  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  std::cout<< "**** Matrix K lidar update ****=  "<<std::endl;
  std::cout<< K <<std::endl;

  //residual
  VectorXd z_diff = z - z_pred;

  //calculate NIS
  NIS_lidar_ = z_diff.transpose() * S.inverse() * z_diff;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;

    //std::cout << "x_after_lidar: " << std::endl << x_ << std::endl;

  P_ = P_ - K*S*K.transpose();
    std::cout<< "**** Matrix P after lidar update ****=  "<<std::endl;
    std::cout<< P_<<std::endl;
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


  // set measurement dimension for radar
  int n_z = 3;

  //create  matrix for sigma point in measurements space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //Transform sigma points into measurement space
  for (int i = 0; i < 2*n_aug_+1; i++) {
    //extract values
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);
    double v   = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    double vx = v * cos(yaw);
    double vy = v * sin(yaw);

    // measurement model
    Zsig(0, i) = sqrt(p_x*p_x + p_y*p_y);
    Zsig(1, i) = atan2(p_y, p_x);
    Zsig(2, i) = (p_x*vx + p_y*vy) / sqrt(p_x*p_x + p_y*p_y);
  }

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);

  for (int i = 0; i < 2*n_aug_ + 1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);

  for (int i = 0; i < 2*n_aug_ + 1; i++){

    VectorXd z_diff(3);
    z_diff = Zsig.col(i) - z_pred;

    //angle normalizatiom
    while (z_diff(1)>M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;
    //z_diff(1) = angle_normalization(z_diff(1));

   // auto addition = weights_(i) * z_diff * z_diff.transpose();
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  MatrixXd R_radar = MatrixXd(3, 3);
  R_radar << std_radr_*std_radr_,                 0,                     0,
          0,             std_radphi_*std_radphi_,                     0,
          0,                                   0, std_radrd_*std_radrd_;

   S = S + R_radar;

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);


  //calculate cross correlation matrix

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1) < M_PI) z_diff(1) += 2.*M_PI;
   // z_diff(1) = angle_normalization(z_diff(1));

    //state differnce
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

      std::cout<< "**** x_diff in radar update ****=  "<<std::endl;
      std::cout<<i<<"  "<< x_diff(3) <<" Xpred= "<< Xsig_pred_.col(i)(3)<<"  x_ = "<<x_(3)<< std::endl;

    //angle normalization
    while (x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;
    //x_diff(3) = angle_normalization(x_diff(3));
      std::cout<<i<<" after angle normalization =  "<< x_diff(3)<<std::endl;


      Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  std::cout<< "**** Matrix Tc radar update ****=  "<<std::endl;
  std::cout<< Tc <<std::endl;


  //Kalman gain K;
   MatrixXd K = Tc * S.inverse();
    std::cout<< "**** Matrix K radar update ****=  "<<std::endl;
    std::cout<< K <<std::endl;

  //extract measurement
  VectorXd z = meas_package.raw_measurements_;

  //residual
  VectorXd z_diff = z - z_pred;


  //update state mean and covariance matrix

  x_ = x_ + K * z_diff;
   // std::cout << "x_after_radar: " << std::endl << x_ << std::endl;
    P_ = P_ - K*S*K.transpose();
    std::cout<< "**** Matrix P after radar ****=  "<<std::endl;
    std::cout<< P_<<std::endl;

  //calculate NIS
  NIS_lidar_ = z_diff.transpose() * S.inverse() * z_diff;

}
