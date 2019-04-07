#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  
  /**
   * TODO: Calculate the RMSE here.
   */
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;
  	
  if (estimations.size() != ground_truth.size() || estimations.size() == 0){
   	std::cout<<"InValid"<<std::endl;
   	return rmse;
  }
    
  for(unsigned int i=0; i<estimations.size(); ++i){
  	VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array()*residual.array();
    rmse+=residual;
  }
  	
  rmse = rmse/estimations.size();
 
  rmse = rmse.array().sqrt();
  
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
  MatrixXd Hj(3, 4);
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);
  
  double d1 = px*px + py*py;
  double d2 = sqrt(d1);
  double d3 = (d1*d2);
  
  double r31 = py*(vx*py - vy*px);
  double r32 = px*(vy*px - vx*py);
  
  // check division by zero
  if (fabs(d1) < 0.0001){
  	std::cout<<"Error - Divide by zero"<<std::endl;
	return Hj;
  }
  // compute the Jacobian matrix
  Hj << px/d2, py/d2, 0, 0,
    	-py/d1, px/d1, 0, 0,
	    r31/d3, r32/d3, px/d2, py/d2;
  return Hj;
}
