#ifndef MEASUREMENT_PACKAGE_H_
#define MEASUREMENT_PACKAGE_H_

#include "Eigen/Dense"
#include <iostream>  // ostream

class MeasurementPackage
{

public:

  // Data members ---------------------------------------------

  long long timestamp_;

  enum SensorType{
    LASER,
    RADAR
  } sensor_type_;

  Eigen::VectorXd raw_measurements_;


  // Member function -------------------------------------------

  // return true if at least one of the raw measurmenets is non-zero
  bool isValid() const;

};


// display MeasurementPackage info
std::ostream& operator<<( std::ostream& os, const MeasurementPackage &meas );




#endif /* MEASUREMENT_PACKAGE_H_ */
