#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF( float noise_ax, float noise_ay )

	:  ekf_(), is_initialized_(false), previous_timestamp_(0),
	   R_laser_( MatrixXd(2, 2) ), R_radar_( MatrixXd(3, 3) ),
	   H_laser_( MatrixXd(2, 4) ), Hj_( MatrixXd(3, 4) ),
	   noise_ax_(noise_ax), noise_ay_(noise_ay)
{

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */

  // measurement covariance matrix - laser (x, y)
  R_laser_ <<  0.0225,  0,
               0,       0.0225;

  // measurement covariance matrix - radar (ro, phi, rod)
  R_radar_ <<  0.09,  0,       0,
               0,     0.0009,  0,
               0,     0,       0.09;

  // measurement matrix - laser
  H_laser_ << 1,  0,  0,  0,
              0,  1,  0,  0;

  // Hj_ is calculated by ProcessMeasurement() at each processing cycle

  // initialize state transition matrix
  // (dt-dependent elements are later updated at each subsequent cycle)
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
			       0, 1, 0, 1,
			       0, 0, 1, 0,
			       0, 0, 0, 1;

	// initialize the process covariance matrix Q
	ekf_.Q_ = MatrixXd(4, 4);

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}






void FusionEKF::ProcessMeasurement( const MeasurementPackage &measurement_pack )
{

	//cout << "Processing measurement...\n" << measurement_pack << endl;

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if ( !is_initialized_ )
  {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */

    // initialize state and state covariance matrix
		ekf_.x_ = VectorXd(4);
		ekf_.x_ << 0, 0, 0, 0;

    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1, 0, 0,    0,
			         0, 1, 0,    0,
			         0, 0, 1000, 0,
			         0, 0, 0, 1000;


  	// case invalid measurement: skip cycle
    if( ! measurement_pack.isValid() )
    {
    	cout << "\nInvalid input measurement, hence initialization skipped" << endl;
    	return;
    }

    // case valid radar measurement
    else if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
    {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      Note that radar cannot observe phi_dot, hence we assume phi_dot==0.
      */
      float ro = measurement_pack.raw_measurements_(0);
      float phi = measurement_pack.raw_measurements_(1);
      float rod = measurement_pack.raw_measurements_(2);

      ekf_.x_(0) = ro*cos(phi);   // x
      ekf_.x_(1) = ro*sin(phi);   // y
      ekf_.x_(2) = rod*cos(phi);  // vx (assuming phi_dot ==0)
      ekf_.x_(3) = rod*sin(phi);  // vy (assuming phi_dot ==0)

      // since radar provides info on radial velocity, reduce initial uncertainty on vx, vy somewhat
      ekf_.P_(2,2) = 100;
      ekf_.P_(3,3) = 100;

    }
    
    // case valid laser measurement
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
    {
      /**
      Initialize state.
      Laser cannot observe velocity, hence we assume vx = vy = 0
      */
      float x = measurement_pack.raw_measurements_(0);
      float y = measurement_pack.raw_measurements_(1);

      ekf_.x_(0) = x;   // x
      ekf_.x_(1) = y;   // y
      ekf_.x_(2) = 0.;  // vx
      ekf_.x_(3) = 0.;  // vy
    }

    // update previous timestamp for next cycle
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;

  	cout << "\nState after initialization:" << endl;
  	cout << "  * x_ = " << ekf_.x_ << endl;
  	cout << "  * P_ = " << ekf_.P_ << endl;

    return;
  }


  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  // if this timestamp is larger than previous timestamp, run prediction step
  // else, skip prediction step
  if( measurement_pack.timestamp_ > previous_timestamp_ )
  {

	  // compute the time elapsed between the current and previous measurements [seconds]
	  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
	  
	  // update previous timestamp for next cycle
	  previous_timestamp_ = measurement_pack.timestamp_;
  
	  // auxiliary variables
	  float dt2 = dt*dt;
	  float dt3 = dt2*dt;
	  float dt4 = dt3*dt;
  
	  // update state-transition matrix with elapsed time
	  ekf_.F_(0, 2) = dt;
	  ekf_.F_(1, 3) = dt;
  
	  // update the process covariance matrix with elapsed time
	  ekf_.Q_ = MatrixXd(4, 4);
	  ekf_.Q_ <<  dt4/4*noise_ax_,  0,               dt3/2*noise_ax_, 0,
	  		        0,                dt4/4*noise_ay_, 0,               dt3/2*noise_ay_,
	  		        dt3/2*noise_ax_,  0,               dt2*noise_ax_,   0,
	  		        0,                dt3/2*noise_ay_, 0,               dt2*noise_ay_;
  
    // run KF prediction step
    ekf_.Predict();
  
    cout << "\nState after prediction:" << endl;
    cout << "  * x_ = " << ekf_.x_ << endl;
    cout << "  * P_ = " << ekf_.P_ << endl;
  }


  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  // if invalid measurement, skip update step
  if( ! measurement_pack.isValid() )
  {
  	cout << "\nInvalid input measurement, hence update step skipped" << endl;
  	return;
  }

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
  {
    // Radar updates

  	// update jacobian of measurement function
    Hj_ = tools.CalculateJacobian( ekf_.x_ );

    // set measurment and measurement covariance matrices to radar values 
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;

    // run EKF measurement update step
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

 	 	cout << "\nState after update from radar:" << endl;
 	 	cout << "  * x_ = " << ekf_.x_ << endl;
 	 	cout << "  * P_ = " << ekf_.P_ << endl;

  } else {
    // Laser updates

    // set measurment and measurement covariance matrices to laser values 
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;

    // run KF measurement update step
    ekf_.Update(measurement_pack.raw_measurements_);

 	 	cout << "\nState after update from laser:" << endl;
 	 	cout << "  * x_ = " << ekf_.x_ << endl;
 	 	cout << "  * P_ = " << ekf_.P_ << endl;

  }

  // print the output
  //cout << "  * x_ = " << ekf_.x_ << endl;
  //cout << "  * P_ = " << ekf_.P_ << endl;
}
