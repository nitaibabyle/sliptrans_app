#include "sliptrans.hpp"

Eigen::Matrix<double, 6, 1> Joint_right; // joint
Eigen::Matrix<double, 6, 1> Joint_left;  // joint
Eigen::Matrix<double, 3, 1> Slip_right;  // slip
Eigen::Matrix<double, 3, 1> Slip_left;   // slip
Eigen::Matrix<double, 3, 1> T;           // torso