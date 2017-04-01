#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;


Tools::Tools() {}

Tools::~Tools() {}


VectorXd Tools::CalculateRMSE(
  const vector<VectorXd> &estimations,
  const vector<VectorXd> &ground_truth)
{
  /**
  TODO:
    * Calculate the RMSE here.
  */
  if( estimations.size() != ground_truth.size() || estimations.size()==0 )
  {
		std::cerr << "Tools::CalculateRMSE: invalid input dimensions.\n" << std::endl;
		exit(EXIT_FAILURE);  	
  }

  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

	//accumulate squared residuals
	for(std::size_t i=0; i<estimations.size(); ++i)
	{
		VectorXd residual = estimations[i] - ground_truth[i];

		//coefficient-wise multiplication
		residual = residual.array() * residual.array();
		rmse += residual;
	}

  return ( rmse/estimations.size() ).array().sqrt();
}


MatrixXd Tools::CalculateJacobian( const VectorXd& x_state ) {

	MatrixXd Hj(3,4);

	//recover state parameters
	float x  = x_state(0);
	float y  = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//pre-compute a set of terms to avoid repeated calculation
	float ro2 = x*x+y*y;
	float ro = std::sqrt(ro2);
	float ro3 = (ro*ro2);

	//check division by zero
	if( ro < 1e-8 ){
		std::cerr << "Tools::CalculateJacobian: division by zero.\n" << std::endl;
		exit(EXIT_FAILURE);
	}

	//compute the Jacobian matrix
	Hj <<  (x/ro),              (y/ro),                0,     0,
		    -(y/ro2),             (x/ro2),               0,     0,
		    y*(vx*y - vy*x)/ro3,  x*(x*vy - y*vx)/ro3,  x/ro,  y/ro;

	return Hj;
}
