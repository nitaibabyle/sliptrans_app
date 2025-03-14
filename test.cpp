#include <iostream>
#include "test.hpp"

fn_slip_trans test1;

fn_slip_trans test2;

int main(int, char **)
{

    T << 5 * 2 * M_PI / 180, 0.0, 0; // rotation matrix

    std::cout << "Hello, from testslip!\n";
    Eigen::Matrix<double, 12, 1> J;  // joint angle
    Eigen::Matrix<double, 12, 1> JV; // joint velocity
    Joint_right << 0, 0, 40 * 2 * M_PI / 180, -80 * 2 * M_PI / 180, 40 * 2 * M_PI / 180, 0;
    Joint_left << 0, 0, 40 * 2 * M_PI / 180, -80 * 2 * M_PI / 180, 40 * 2 * M_PI / 180, 0;
    J.block<6, 1>(0, 0) = Joint_right;
    J.block<6, 1>(6, 0) = Joint_left;
    JV.block<6, 1>(0, 0) = 0 * Joint_right;
    JV.block<6, 1>(6, 0) = 0 * Joint_left;
    for (int i = 0; i < 1000; i++)
    {
        test1.jointAngel_input = J;
        test1.jointAngelV_input = JV;
        test1.torso_input = T;
        test1.calculate_slipstate();
        test1.slipState_output;
        test1.slipStateV_output;
        std::cout << test1.slipState_output << std::endl;

        test2.slipState_input = test1.slipState_output;
        test2.slipStateV_input = test1.slipStateV_output;
        test2.torso_input = T;
        test2.calculate_jointangle();
        test2.jointAngel_output;
        test2.jointAngelV_output;
        std::cout << test2.jointAngel_output * 90 / M_PI << std::endl;
    }

    return 0;
}
