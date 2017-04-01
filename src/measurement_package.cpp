#include "measurement_package.h"

#include <iostream>  // std::ostream

using Eigen::VectorXd;


std::ostream& operator<<( std::ostream& os, const MeasurementPackage &meas )  
{  
  os << "\nMeasurementPackage" << std::endl;
  os << "  * time stamp : " << meas.timestamp_ << std::endl;
  os << "  * sensor type : " << ( meas.sensor_type_==MeasurementPackage::LASER ? "laser" : "radar" ) << std::endl;
  os << "  * raw measurements : " << meas.raw_measurements_ << std::endl;
  return os;  
}


bool MeasurementPackage::isValid() const
{
	for( size_t i=0; i<raw_measurements_.size(); ++i )
	{
		if( raw_measurements_(i) != 0 )
			return true;
	}

	return false;
}
