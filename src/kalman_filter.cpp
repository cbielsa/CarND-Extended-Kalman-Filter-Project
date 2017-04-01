#include "kalman_filter.h"
#include <cmath>     // std::sqrt()
#include <iostream>  // std::cerr

using Eigen::MatrixXd;
using Eigen::VectorXd;


KalmanFilter::KalmanFilter()
	: I_( MatrixXd::Identity(4,4) )
 {}


KalmanFilter::~KalmanFilter() {}


void KalmanFilter::Init(
	const VectorXd &x_in, const MatrixXd &P_in, const MatrixXd &F_in,
    const MatrixXd &H_in, const MatrixXd &R_in, const MatrixXd &Q_in )
{
    x_ = x_in;
    P_ = P_in;
    F_ = F_in;
    H_ = H_in;
	R_ = R_in;
	Q_ = Q_in;
}


void KalmanFilter::Predict()
{
	x_ = F_ * x_;
	P_ = F_ * P_ * F_.transpose() + Q_;
	return;
}


void KalmanFilter::Update( const VectorXd &z )
{
	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd K = P_ * Ht * S.inverse();
	
	// new estimate
	x_ += (K * y);
	P_ = (I_ - K * H_) * P_;
	return;
}


VectorXd KalmanFilter::MeasurementFunctionRadar() const
{
	VectorXd vZexp(3);

	float x  = x_(0);
	float y  = x_(1);
	float vx = x_(2);
	float vy = x_(3);

	float ro = std::sqrt(x*x + y*y);
	if( ro < 1e-8 )
	{
		std::cerr << "KalmanFilter::MeasurementFunctionRadar: division by zero.\n" << std::endl;
		exit(EXIT_FAILURE);
	}

	vZexp << ro, std::atan2(y,x), (x*vx+y*vy)/ro;
	return vZexp;
}


void KalmanFilter::UpdateEKF( const VectorXd &z )
{

	// calculate measurements consistent with current state estimate
	VectorXd z_pred = MeasurementFunctionRadar();

	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd K = P_ * Ht * S.inverse();
	
	// new estimate
	x_ += (K * y);
	P_ = (I_ - K * H_) * P_;
	return;
}



