#include <iostream>
#include "hlqIkQPbased.hpp"
#include <iomanip>
#include <iostream>

/*
	Eigen documation https://eigen.tuxfamily.org/dox/group__TutorialMatrixClass.html
*/

ikQP::ikQP()
{

	zzx0.resize(8); // right4 left4
	zzx0.setZero();

	zzxResult.resize(8 * 2); // angle8 velocity8
	zzxResult.setZero();

	zzW1.resize(3); // right weight
	zzW1.diagonal() << 10, 10, 10;

	zzW2.resize(3); // legt weight
	zzW2.diagonal() << 10, 10, 10;

	zzW3.resize(8);
	// expect qdot smaller
	zzW3.diagonal() << 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01;

	zzW4.resize(8); // expect qdot changes smaller
	zzW4.diagonal() << 1, 1, 1, 1, 1, 1, 1, 1;

	zzW4.resize(8); // yaw velocity
	zzW4.diagonal() << 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01;

	zzW5.resize(2);
	zzW5.diagonal() << 200, 200;

	zzM1.resize(3, 8);
	zzM1.setZero();
	zzb1.resize(3, 1);
	zzb1.setZero();
	zzM2.resize(3, 8);
	zzM2.setZero();
	zzb2.resize(3, 1);
	zzb2.setZero();
	zzM3.resize(8, 8);
	zzM3 << Eigen::MatrixXd::Identity(8, 8);
	zzb3.resize(8, 1);
	zzb3.setZero();
	zzM4.resize(8, 8);
	zzM4 << Eigen::MatrixXd::Identity(8, 8);
	zzb4.resize(8, 1);
	zzb4.setZero();
	zzM5.resize(2, 8);
	zzM5 << 0, 1, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 1, 0, 0;
	zzb5.resize(2, 1);
	zzb5.setZero(); // yaw set to zero, later may change and control

	// zzM5.resize(2, 6);
	// zzM5 << 0, 1, 1, 1, 0, 0, 0, 0,
	// 	0, 0, 0, 0, 0, 1, 1, 1;
	// zzb5.resize(2, 1);
	// zzb5.setZero();

	zzAineq1.resize(14, 8);
	zzAineq1 << 0, 0, 0, 1, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 1,
		1, 0, 0, 0, 0, 0, 0, 0,
		-1, 0, 0, 0, 0, 0, 0, 0,
		0, 1, 0, 0, 0, 0, 0, 0,
		0, -1, 0, 0, 0, 0, 0, 0,
		0, 0, 1, 0, 0, 0, 0, 0,
		0, 0, -1, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 1, 0, 0, 0,
		0, 0, 0, 0, -1, 0, 0, 0,
		0, 0, 0, 0, 0, 1, 0, 0,
		0, 0, 0, 0, 0, -1, 0, 0,
		0, 0, 0, 0, 0, 0, 1, 0,
		0, 0, 0, 0, 0, 0, -1, 0;

	zzbineq1.resize(14);
	zzbineq1.setZero();

	zzAineq2.resize(8 * 2, 8);
	zzAineq2.block(0, 0, 8, 8) = Eigen::MatrixXd::Identity(8, 8);
	zzAineq2.block(8, 0, 8, 8) = -Eigen::MatrixXd::Identity(8, 8);
	zzbineq2.resize(8 * 2);
	zzbineq2 << (DEGREE_2_RAD * ANGLE_STEP_MAX) * Eigen::VectorXd::Ones(8 * 2);

	zzAeq1.resize(2, 8);
	zzAeq1 << 0, 1, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 1, 0, 0;

	zzbeq1.resize(2);
	zzbeq1.setZero();

	// 初始化QP问题
	zzqpSolver.setObjectFuncParam<hlqGeneralQPForm::mo111>(zzM1, zzb1, zzW1);
	zzqpSolver.setObjectFuncParam<hlqGeneralQPForm::mo111>(zzM2, zzb2, zzW2);
	zzqpSolver.setObjectFuncParam<hlqGeneralQPForm::mo111>(zzM3, zzb3, zzW3);
	zzqpSolver.setObjectFuncParam<hlqGeneralQPForm::mo111>(zzM4, zzb4, zzW4);
	zzqpSolver.setObjectFuncParam<hlqGeneralQPForm::mo111>(zzM5, zzb5, zzW5);
	zzqpSolver.setInequationConstraint<hlqGeneralQPForm::mo111>(zzAineq1, zzbineq1);
	zzqpSolver.setInequationConstraint<hlqGeneralQPForm::mo111>(zzAineq2, zzbineq2);
	zzqpSolver.setEquationConstraint<hlqGeneralQPForm::mo111>(zzAeq1, zzbeq1);
	zzqpSolver.clearNumTimer();
}

ikQP::~ikQP()
{
}

void ikQP::slipgetikResult(robot_para &ref_robot, Eigen::VectorXd &targetPosIn, Eigen::VectorXd &targetVelIn, Eigen::VectorXd &targetAngleIn, Eigen::VectorXd &targetYawvIn, double stepTimeIn)
{
	// converge coeffiecient
	double dKp = 10.0;
	// set local qref current
	RLeg_jointAngle = ref_robot.jointAngle.head(6);
	LLeg_jointAngle = ref_robot.jointAngle.tail(6);
	RLeg_jointAngleV = ref_robot.jointAngleV.head(6);
	LLeg_jointAngleV = ref_robot.jointAngleV.tail(6);
	zzx0 << RLeg_jointAngle.head(4), LLeg_jointAngle.head(4);

	// qref changes is the opt valuables, or qdot*T

	// target, the difference between the target and current qref
	Eigen::Vector3d rtargetVel = targetVelIn.head(3) + dKp * (targetPosIn.head(3) - ref_robot.contactPointPos_right.block(0, 0, 3, 1));
	Eigen::Vector3d ltargetVel = targetVelIn.tail(3) + dKp * (targetPosIn.tail(3) - ref_robot.contactPointPos_left.block(0, 0, 3, 1));

	// right foot end effector target 3
	zzM1.block(0, 0, 3, 4) = ref_robot.linkJacobian_right.block<3, 4>(0, 0) / stepTimeIn;
	zzb1 << rtargetVel;

	zzqpSolver.setObjectFuncParam<hlqGeneralQPForm::mo110>(zzM1, zzb1, zzW1);

	// //legt foot end effector target 3
	zzM2.block(0, 4, 3, 4) = ref_robot.linkJacobian_left.block<3, 4>(0, 0) / stepTimeIn;
	zzb2 << ltargetVel;
	zzqpSolver.setObjectFuncParam<hlqGeneralQPForm::mo110>(zzM2, zzb2, zzW2);

	// qref's changes,or qdot need to small enough 8
	zzqpSolver.setObjectFuncParam<hlqGeneralQPForm::mo000>(zzM3, zzb3, zzW3);

	// qref's changing need to similar enough with last time  8
	zzb4 << RLeg_jointAngleV.head(4) * stepTimeIn, LLeg_jointAngleV.head(4) * stepTimeIn;
	zzqpSolver.setObjectFuncParam<hlqGeneralQPForm::mo010>(zzM4, zzb4, zzW4);

	// yaw small 2
	// zzb5 << (targetYawvIn[0] - zzx0[1]) / stepTimeIn, (targetYawvIn[1] - zzx0[5]) / stepTimeIn;
	// zzb5 << targetYawvIn.head(2);
	// zzqpSolver.setObjectFuncParam<hlqGeneralQPForm::mo010>(zzM5, zzb5, zzW5);

	// 先让脚板与地面尽量保持平衄1�7 2
	// 22.02.14修改为可变角庄1�7 //zjt not consider ankle
	// b5 << -x0[1] - x0[2] - x0[3] + targetAngleIn[0],
	// 	-x0[5] - x0[6] - x0[7] + targetAngleIn[1];
	// zzqpSolver.setObjectFuncParam<hlqGeneralQPForm::mo010>(M5, b5, W5);

	/*  等式约束函数  */

	/*  不等式约束函敄1�7  */
	// knee negative constrant, qref limited           14
	zzbineq1 << 0 - zzx0[3], // right knee smaller than 0
		0 - zzx0[7],		 // left knee
		DEGREE_2_RAD * (ANKLE_ANGLE_MAX / 1) - zzx0[0],
		DEGREE_2_RAD * (ANKLE_ANGLE_MAX / 1) + zzx0[0],
		DEGREE_2_RAD * (ANKLE_ANGLE_MAX / 4) - zzx0[1],
		DEGREE_2_RAD * (ANKLE_ANGLE_MAX / 4) + zzx0[1],
		DEGREE_2_RAD * (ANKLE_ANGLE_MAX)-zzx0[2],
		DEGREE_2_RAD * (ANKLE_ANGLE_MAX) + zzx0[2],
		DEGREE_2_RAD * (ANKLE_ANGLE_MAX / 1) - zzx0[4],
		DEGREE_2_RAD * (ANKLE_ANGLE_MAX / 1) + zzx0[4],
		DEGREE_2_RAD * (ANKLE_ANGLE_MAX / 4) - zzx0[5],
		DEGREE_2_RAD * (ANKLE_ANGLE_MAX / 4) + zzx0[5],
		DEGREE_2_RAD * (ANKLE_ANGLE_MAX)-zzx0[6],
		DEGREE_2_RAD * (ANKLE_ANGLE_MAX) + zzx0[6];
	zzqpSolver.setInequationConstraint<hlqGeneralQPForm::mo010>(zzAineq1, zzbineq1);

	// qdot*T limited   16
	zzqpSolver.setInequationConstraint<hlqGeneralQPForm::mo000>(zzAineq2, zzbineq2);

	// hip yaw constrain to zero   2
	zzqpSolver.setEquationConstraint<hlqGeneralQPForm::mo000>(zzAeq1, zzbeq1);

	// calculate result qref changes
	zzqpSolver.quadprogSolver();

	// integrate the result of qref changes to gain new qref
	zzxResult.head(8) = zzx0 + zzqpSolver.qpState;

	// std::cout << "\n"
	// 		  << zzqpSolver.qpState << "\n"
	// 		  << std::endl;

	// new qdot
	zzxResult.tail(8) = zzqpSolver.qpState / stepTimeIn;

	// set ref traj
	ref_robot.jointAngle.block<4, 1>(0, 0) << zzxResult[0], zzxResult[1], zzxResult[2], zzxResult[3];
	ref_robot.jointAngleV.block<4, 1>(0, 0) << zzxResult[8], zzxResult[9], zzxResult[10], zzxResult[11];
	ref_robot.jointAngle.block<4, 1>(6, 0) << zzxResult[4], zzxResult[5], zzxResult[6], zzxResult[7];
	ref_robot.jointAngleV.block<4, 1>(6, 0) << zzxResult[12], zzxResult[13], zzxResult[14], zzxResult[15];
}
