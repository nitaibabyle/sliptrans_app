#ifndef HLQ_IK_QPBASED_HPP
#define HLQ_IK_QPBASED_HPP

#include "hlqGeneralQPForm.hpp"
#include <Eigen/Dense>
#include "../robotpara.hpp"

#define ANKLE_ANGLE_MAX 80.0
#define ANGLE_STEP_MAX 1
#define DEGREE_2_RAD 0.03490658

class ikQP
{
public:
	ikQP();
	~ikQP();

	Eigen::VectorXd zzx0;
	Eigen::VectorXd zzxResult;

	// ��ￄ1�7?
	void slipgetikResult(robot_para &robot, Eigen::VectorXd &targetPosIn, Eigen::VectorXd &targetVelIn, Eigen::VectorXd &targetAngleIn, Eigen::VectorXd &targetYawvIn, double stepTimeIn);

private:
	Eigen::Matrix<double, 6, 1> RLeg_jointAngle;
	Eigen::Matrix<double, 6, 1> LLeg_jointAngle;
	Eigen::Matrix<double, 6, 1> RLeg_jointAngleV;
	Eigen::Matrix<double, 6, 1> LLeg_jointAngleV;
	// QP��׼����ￄ1�7?

	//      state number, objective number, equality number, inequality number
	hlqGeneralQPForm::quadprog<int, 8, 24, 4, 30> zzqpSolver;

	Eigen::DiagonalMatrix<double, Eigen::Dynamic> zzW1;
	Eigen::DiagonalMatrix<double, Eigen::Dynamic> zzW2;
	Eigen::DiagonalMatrix<double, Eigen::Dynamic> zzW3;
	Eigen::DiagonalMatrix<double, Eigen::Dynamic> zzW4;
	Eigen::DiagonalMatrix<double, Eigen::Dynamic> zzW5;
	//  目锟疥函锟斤拷锟侥撅拷锟斤拷
	Eigen::MatrixXd zzM1;
	Eigen::VectorXd zzb1;
	Eigen::MatrixXd zzM2;
	Eigen::VectorXd zzb2;
	Eigen::MatrixXd zzM3;
	Eigen::VectorXd zzb3;
	Eigen::MatrixXd zzM4;
	Eigen::VectorXd zzb4;
	Eigen::MatrixXd zzM5;
	Eigen::VectorXd zzb5;
	//  锟斤拷式约锟斤拷锟侥撅拷锟斤拄1�7

	//  锟斤拷锟斤拷式约锟斤拷锟侥撅拷锟斤拷
	Eigen::MatrixXd zzAineq1;
	Eigen::VectorXd zzbineq1;
	Eigen::MatrixXd zzAineq2;
	Eigen::VectorXd zzbineq2;

	Eigen::MatrixXd zzAeq1;
	Eigen::VectorXd zzbeq1;
};

#endif // HLQ_IK_QPBASED_HPP
