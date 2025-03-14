#include <Eigen/Dense>

class robot_para
{
public:
    robot_para();
    ~robot_para();
    //  hiproll hipyaw hippitch kneepitch anklepitch ankleroll
    Eigen::Matrix<double, 12, 1> jointAngle;  // right 6 left 6
    Eigen::Matrix<double, 12, 1> jointAngleV; // right 6 left 6

    //   length  rollangle  pitchangle
    Eigen::Matrix<double, 6, 1> slipState;  // right 3 left 3
    Eigen::Matrix<double, 6, 1> slipStateV; // right 3 left 3

    Eigen::Matrix<double, 3, 1> slipfootPos_right; // from joint space to foot end space, forward kinematic, used in act foot pos calculation
    Eigen::Matrix<double, 3, 1> slipfootVel_right; // from joint space to foot end space, forward kinematic, used in act foot pos calculation
    Eigen::Matrix<double, 3, 1> slipfootPos_left;  // from joint space to foot end space, forward kinematic, used in act foot pos calculation
    Eigen::Matrix<double, 3, 1> slipfootVel_left;  // from joint space to foot end space, forward kinematic, used in act foot pos calculation

    Eigen::Vector2d slip_leg_length;      // from foot end space to slip space. forward kinematic, used in act slip state calculation
    Eigen::Vector2d slip_leg_roll_angle;  // from foot end space to slip space. forward kinematic, used in act slip state calculation
    Eigen::Vector2d slip_leg_pitch_angle; // from foot end space to slip space. forward kinematic, used in act slip state calculation

    Eigen::Matrix<double, 3, 1> foot_pos_slip_right; // from slip space to foot end pos, inverse kinematic, used in reference foot pos calculation
    Eigen::Matrix<double, 3, 1> foot_vel_slip_right; // from slip space to foot end pos, inverse kinematic, used in reference foot pos calculation
    Eigen::Matrix<double, 3, 1> foot_pos_slip_left;  // from slip space to foot end pos, inverse kinematic, used in reference foot pos calculation
    Eigen::Matrix<double, 3, 1> foot_vel_slip_left;  // from slip space to foot end pos, inverse kinematic, used in reference foot pos calculation

    Eigen::Matrix<double, 3, 1> contactPointPos_left;  // from foot end space to joint space. inverse kinematic, used in reference joint angle calculation
    Eigen::Matrix<double, 3, 1> contactPointPos_right; // from foot end space to joint space. inverse kinematic, used in reference joint angle calculation
    Eigen::Matrix<double, 3, 1> contactPointVel_left;  // from foot end space to joint space. inverse kinematic, used in reference joint angle calculation
    Eigen::Matrix<double, 3, 1> contactPointVel_right; // from foot end space to joint space. inverse kinematic, used in reference joint angle calculation

    Eigen::Matrix<double, 3, 6> linkJacobian_right; // ÿ�����˵������ſɱ�
    Eigen::Matrix<double, 3, 6> linkJacobian_left;  // ÿ�����˵������ſɱ�
};