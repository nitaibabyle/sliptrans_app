#ifndef HLQ_GENERAL_QP_FORM_HPP
#define HLQ_GENERAL_QP_FORM_HPP
/*
	min 1/2*||Mx-b||^2

	!BASED
	min 1/2*x^T*H*x+f^T*x
	s.t. Aeq*x = beq
		   A*x < b
	Matrix <double, nrvar, nrvar> H;
	Matrix <double, nrvar,1 > f;
	Matrix <double, nreq, nrvar> Aeq;
	Matrix <double, nreq, 1> Beq;
	Matrix <double, nrineq, nrvar> Aineq;
	Matrix <double, nrineq, 1> Bineq;
*/
#include <iostream>
#include <Eigen/Dense>
#include <memory>
#include "QP_Pierre/QuadProg.h"

#pragma warning(disable : 4244)

namespace hlqGeneralQPForm
{

	/// ���²���ʱ��ģʽѡ��  0 0 0
	const int mo000 = 0; // ��������
	const int mo001 = 1; // ��3λ����
	const int mo010 = 2; // ��2λ
	const int mo011 = 3; // ��2��3λ
	const int mo100 = 4; // ��1λ
	const int mo101 = 5; // ��1��3λ
	const int mo110 = 6; // ��1��2λ
	const int mo111 = 7; // ��1��2��3λ

	template <class T, const int SNum = 0, const int ONum = 0, const int EqNum = 0, const int IneqNum = 0>
	class quadprog
	{
	public:
		quadprog();
		~quadprog();

		// ���Ŀ�꺯��?
		template <const int mode>
		bool setObjectFuncParam(Eigen::MatrixXd &subMIn, Eigen::VectorXd &subbIn, Eigen::DiagonalMatrix<double, Eigen::Dynamic> &subWIn);

		// ֱ������Ŀ�꺯��
		bool setUpdateObjectFuncParam(Eigen::MatrixXd &MIn, Eigen::VectorXd &bIn);

		// ��ϵ�ʽԼ��?
		template <const int mode>
		bool setEquationConstraint(Eigen::MatrixXd &subAeqIn, Eigen::VectorXd &subbeqIn);

		// ��ϲ���ʽԼ��?
		template <const int mode>
		bool setInequationConstraint(Eigen::MatrixXd &subAIn, Eigen::VectorXd &subbIn);

		// �Ż����?
		bool quadprogSolver();

		// ���Ĳ���Ŀ�����?
		bool updateObjectMatrix();

		// ���Ĳ��ֵ�ʽ����
		bool updateEquationMatrix();

		// ���Ĳ��ֲ���ʽ����
		bool updateInequationMatrix();

		// ��ʾ����
		void displayQPMatrix();

		// ���������?
		void clearNumTimer();

		// �Ż���������
		int qpStateNum = SNum;

		// �Ż�����
		Eigen::VectorXd qpState;

		// Ŀ�꺯������ ����������
		int qpObjectFuncNum;

		// ��ʽԼ���������� ����������
		int qpEqConstraintNum;

		// ����ʽԼ���������� ����������
		int qpIneqConstraintNum;

		// �Ƿ�ֱ���趨��H f �����־�? Ĭ��=0 ʹ��M*x-b��һ�׵�
		int qpHfDirectSetFlag;

	private:
		// ʹ������ָ���� C2280 ����������ɾ���ĺ��� �������?
		std::shared_ptr<Eigen::QuadProgDense> solver = std::make_shared<Eigen::QuadProgDense>(SNum, EqNum, IneqNum);

		// ����QP�Ż���ʽϵ������
		Eigen::MatrixXd M;
		// ����QP�Ż���ʽĿ������
		Eigen::VectorXd bm;
		// ����QP�Ż���ʽȨ��ϵ��
		Eigen::DiagonalMatrix<double, Eigen::Dynamic> W;

		// Ŀ�� Hessian����
		Eigen::MatrixXd H;
		// �ݶ�ϵ������
		Eigen::VectorXd f;
		// ��ʽϵ������
		Eigen::MatrixXd Aeq;
		// ��ʽ����
		Eigen::VectorXd beq;

		// ����ʽϵ������
		Eigen::MatrixXd A;
		// ����ʽ����
		Eigen::VectorXd b;
	};

	template <class T, const int SNum, const int ONum, const int EqNum, const int IneqNum>
	quadprog<T, SNum, ONum, EqNum, IneqNum>::quadprog() : qpObjectFuncNum(0), qpEqConstraintNum(0), qpIneqConstraintNum(0), qpHfDirectSetFlag(0)
	{
		qpState.resize(SNum);
		qpState.setZero();

		Aeq.resize(EqNum, SNum);
		Aeq.setZero();
		beq.resize(EqNum);
		beq.setZero();

		A.resize(IneqNum, SNum);
		A.setZero();
		b.resize(IneqNum);
		b.setZero();

		M.resize(ONum, SNum);
		M.setZero();
		bm.resize(ONum);
		bm.setZero();
		W.resize(ONum);
		W.setZero();

		H.resize(SNum, SNum);
		H.setZero();
		f.resize(SNum);
		f.setZero();
	}

	template <class T, const int SNum, const int ONum, const int EqNum, const int IneqNum>
	quadprog<T, SNum, ONum, EqNum, IneqNum>::~quadprog()
	{
	}

	// ���Ŀ�꺯��?
	template <class T, const int SNum, const int ONum, const int EqNum, const int IneqNum>
	template <const int mode>
	bool quadprog<T, SNum, ONum, EqNum, IneqNum>::setObjectFuncParam(Eigen::MatrixXd &subMIn, Eigen::VectorXd &subbIn, Eigen::DiagonalMatrix<double, Eigen::Dynamic> &subWIn)
	{
		// ��øô����ӵľ����С
		int nObjectNum = subMIn.rows();

		if ((this->qpObjectFuncNum + nObjectNum) > ONum)
		{
			std::cout << "The object number is wrong!" << std::endl;
			return 1;
		}

		switch (mode)
		{
		case mo000:
			break;

		case mo010:
		{
			this->bm.block(this->qpObjectFuncNum, 0, nObjectNum, 1) = subWIn * subbIn;
			break;
		}
		case mo100:
		{
			this->M.block(this->qpObjectFuncNum, 0, nObjectNum, SNum) = subWIn * subMIn;
			break;
		}
		case mo001:
		case mo011:
		case mo101:
		case mo110:
		case mo111:
		{
			this->M.block(this->qpObjectFuncNum, 0, nObjectNum, SNum) = subWIn * subMIn;
			this->bm.block(this->qpObjectFuncNum, 0, nObjectNum, 1) = subWIn * subbIn;
			break;
		}
		default:
			break;
		}

		this->qpObjectFuncNum += nObjectNum;

		return 0;
	}

	// ���Ŀ�꺯��?
	template <class T, const int SNum, const int ONum, const int EqNum, const int IneqNum>
	bool quadprog<T, SNum, ONum, EqNum, IneqNum>::setUpdateObjectFuncParam(Eigen::MatrixXd &MIn, Eigen::VectorXd &bIn)
	{

		this->H = MIn;
		this->f = bIn;

		this->qpObjectFuncNum = ONum - 1;

		this->qpHfDirectSetFlag = 1; // ��ʾ����û����M*x-bm��һ�� ���ü���H f

		return 0;
	}

	// �������Լ��?
	template <class T, const int SNum, const int ONum, const int EqNum, const int IneqNum>
	template <const int mode>
	bool quadprog<T, SNum, ONum, EqNum, IneqNum>::setEquationConstraint(Eigen::MatrixXd &subAeqIn, Eigen::VectorXd &subbeqIn)
	{

		// ��øô����ӵľ����С
		int nConstraintNum = subAeqIn.rows();

		if ((this->qpEqConstraintNum + nConstraintNum) > EqNum)
		{
			std::cout << "The equation constrain number is wrong!" << std::endl;
			return 1;
		}

		switch (mode)
		{
		case mo000:
		case mo001:
			break;
		case mo010:
		case mo011:
			this->beq.block(this->qpEqConstraintNum, 0, nConstraintNum, 1) = subbeqIn;
			break;
		case mo100:
		case mo101:
			this->Aeq.block(this->qpEqConstraintNum, 0, nConstraintNum, SNum) = subAeqIn;
			break;
		case mo110:
		case mo111:
			this->beq.block(this->qpEqConstraintNum, 0, nConstraintNum, 1) = subbeqIn;
			this->Aeq.block(this->qpEqConstraintNum, 0, nConstraintNum, SNum) = subAeqIn;
			break;
		default:
			break;
		}

		this->qpEqConstraintNum += nConstraintNum;

		return 0;
	}

	template <class T, const int SNum, const int ONum, const int EqNum, const int IneqNum>
	bool quadprog<T, SNum, ONum, EqNum, IneqNum>::updateObjectMatrix()
	{

		return 0;
	}

	template <class T, const int SNum, const int ONum, const int EqNum, const int IneqNum>
	bool quadprog<T, SNum, ONum, EqNum, IneqNum>::updateEquationMatrix()
	{

		return 0;
	}

	template <class T, const int SNum, const int ONum, const int EqNum, const int IneqNum>
	bool quadprog<T, SNum, ONum, EqNum, IneqNum>::updateInequationMatrix()
	{

		return 0;
	}

	// �������Լ��?
	template <class T, const int SNum, const int ONum, const int EqNum, const int IneqNum>
	template <const int mode>
	bool quadprog<T, SNum, ONum, EqNum, IneqNum>::setInequationConstraint(Eigen::MatrixXd &subAIn, Eigen::VectorXd &subbIn)
	{
		// ��øô����ӵľ����С
		int nConstraintNum = subAIn.rows();

		if ((this->qpIneqConstraintNum + nConstraintNum) > IneqNum)
		{
			std::cout << "The inequation constrain number is wrong!" << std::endl;
			return 1;
		}

		switch (mode)
		{
		case mo000:
		case mo001:
			break;
		case mo010:
		case mo011:
			this->b.block(this->qpIneqConstraintNum, 0, nConstraintNum, 1) = subbIn;
			break;
		case mo100:
		case mo101:
			this->A.block(this->qpIneqConstraintNum, 0, nConstraintNum, SNum) = subAIn;
			break;
		case mo110:
		case mo111:
			this->A.block(this->qpIneqConstraintNum, 0, nConstraintNum, SNum) = subAIn;
			this->b.block(this->qpIneqConstraintNum, 0, nConstraintNum, 1) = subbIn;
			break;
		default:
			break;
		}

		this->qpIneqConstraintNum += nConstraintNum;

		return 0;
	}

	// �Ż����?
	template <class T, const int SNum, const int ONum, const int EqNum, const int IneqNum>
	bool quadprog<T, SNum, ONum, EqNum, IneqNum>::quadprogSolver()
	{

		if (this->qpHfDirectSetFlag == 0)
		{

			this->H.triangularView<Eigen::Upper>() = this->M.transpose() * this->M; //
			this->H.triangularView<Eigen::Lower>() = this->H.triangularView<Eigen::Upper>().transpose();

			this->f = -this->M.transpose() * this->bm; //
		}
		else
		{
			this->qpHfDirectSetFlag = 0;
		}

		// solve the QP problem
		int solverResult = solver->solve(this->H, this->f, this->Aeq, this->beq, this->A, this->b);

		// get the controller input
		this->qpState = solver->result();

		// clear
		this->qpObjectFuncNum = 0;
		this->qpEqConstraintNum = 0;
		this->qpIneqConstraintNum = 0;

		return solverResult;
	}

	template <class T, const int SNum, const int ONum, const int EqNum, const int IneqNum>
	void quadprog<T, SNum, ONum, EqNum, IneqNum>::clearNumTimer()
	{
		// clear
		this->qpObjectFuncNum = 0;
		this->qpEqConstraintNum = 0;
		this->qpIneqConstraintNum = 0;
	}

	template <class T, const int SNum, const int ONum, const int EqNum, const int IneqNum>
	void quadprog<T, SNum, ONum, EqNum, IneqNum>::displayQPMatrix()
	{
		std::cout << "********************************************" << std::endl;
		std::cout << "QP form H: \n"
				  << this->H << std::endl;
		std::cout << "QP form f: \n"
				  << this->f << std::endl;
		std::cout << "QP form Aeq: \n"
				  << this->Aeq << std::endl;
		std::cout << "QP form beq: \n"
				  << this->beq << std::endl;
		std::cout << "QP form A: \n"
				  << this->A << std::endl;
		std::cout << "QP form b: \n"
				  << this->b << std::endl;
		std::cout << "********************************************" << std::endl;
	}
}

#endif // !HLQ_GENERAL_QP_FORM_HPP
