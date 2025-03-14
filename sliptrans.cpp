
#include "sliptrans.hpp"

using namespace pinocchio;

fn_slip_trans::fn_slip_trans()
{
    joint_angle.resize(14, 1);
    joint_angleIK.resize(14, 1);

    joint_vel.resize(14, 1);
    joint_velIK.resize(14, 1);

    J_r.resize(6, 20);
    J_l.resize(6, 20);
    // J_.resize(6, 16);

    pinocchio::urdf::buildModel(urdf_file, model_);
    pinocchio::urdf::buildModel(urdf_file_IK, model_IK);

    data_ = pinocchio::Data(model_);
    data_IK = pinocchio::Data(model_IK);

    q_ = pinocchio::neutral(model_);
    q_IK = pinocchio::neutral(model_IK);

    q_ikref = pinocchio::neutral(model_IK);
    qd_ = pinocchio::neutral(model_).block<20, 1>(0, 0);
    qd_IK = pinocchio::neutral(model_IK).block<20, 1>(0, 0);

    q_ikref.setZero();

    Eigen::Matrix<double, 7, 1> floatingBase;
    floatingBase << 0, 0, 0, 0, 0, 0, 1;
    q_.block<7, 1>(0, 0) = floatingBase;
    qd_.block<6, 1>(0, 0) = floatingBase.block<6, 1>(0, 0);

    q_IK.block<7, 1>(0, 0) = floatingBase;
    qd_IK.block<6, 1>(0, 0) = floatingBase.block<6, 1>(0, 0);

    joint_angle.setZero();
    joint_vel.setZero();
    Lfoot_id = model_.getFrameId("Lfoot");
    Rfoot_id = model_.getFrameId("Rfoot");
    torso_id = model_.getFrameId("torso");

    Lfoot_id_IK = model_IK.getFrameId("Lfoot");
    Rfoot_id_IK = model_IK.getFrameId("Rfoot");
    torso_id_IK = model_IK.getFrameId("torso");

    rfoot_joint_id = model_.getJointId("rfootpitch");
    lfoot_joint_id = model_.getJointId("lfootpitch");

    rfoot_joint_id_IK = model_IK.getJointId("rfootpitch");
    lfoot_joint_id_IK = model_IK.getJointId("lfootpitch");

    J_r.setZero();
    J_l.setZero();
    // J_.setZero();

    std::cout << "model name: " << model_.name << std::endl;
    std::cout << "joint number: " << model_.nq << std::endl;
}

fn_slip_trans::~fn_slip_trans()
{
}

void fn_slip_trans::robot2slip(robot_para &act, Eigen::Matrix<double, 3, 1> &torsoAngle)
{
    double leg_hip_roll_angle_right = act.jointAngle(0);
    double leg_hip_yaw_angle_right = act.jointAngle(1);
    double leg_hip_pitch_angle_right = act.jointAngle(2);
    double leg_knee_angle_right = act.jointAngle(3);
    double dleg_hip_roll_angle_right = act.jointAngleV(0);
    double dleg_hip_yaw_angle_right = act.jointAngleV(1);
    double dleg_hip_pitch_angle_right = act.jointAngleV(2);
    double dleg_knee_angle_right = act.jointAngleV(3);

    double leg_hip_roll_angle_left = act.jointAngle(6);
    double leg_hip_yaw_angle_left = act.jointAngle(7);
    double leg_hip_pitch_angle_left = act.jointAngle(8);
    double leg_knee_angle_left = act.jointAngle(9);
    double dleg_hip_roll_angle_left = act.jointAngleV(6);
    double dleg_hip_yaw_angle_left = act.jointAngleV(7);
    double dleg_hip_pitch_angle_left = act.jointAngleV(8);
    double dleg_knee_angle_left = act.jointAngleV(9);

    /*set the kinematic value to pinocchio*/
    joint_angle(0) = leg_hip_roll_angle_left;
    joint_angle(1) = leg_hip_yaw_angle_left;
    joint_angle(2) = leg_hip_pitch_angle_left;
    joint_angle(3) = leg_knee_angle_left;
    joint_angle(4) = 0;
    joint_angle(5) = 0;
    joint_angle(6) = 0;
    joint_angle(7) = leg_hip_roll_angle_right;
    joint_angle(8) = leg_hip_yaw_angle_right;
    joint_angle(9) = leg_hip_pitch_angle_right;
    joint_angle(10) = leg_knee_angle_right;
    joint_angle(11) = 0;
    joint_angle(12) = 0;
    joint_angle(13) = 0;

    joint_vel(0) = dleg_hip_roll_angle_left;
    joint_vel(1) = dleg_hip_yaw_angle_left;
    joint_vel(2) = dleg_hip_pitch_angle_left;
    joint_vel(3) = dleg_knee_angle_left;
    joint_vel(4) = 0.0;
    joint_vel(5) = 0.0;
    joint_vel(6) = 0.0;
    joint_vel(7) = dleg_hip_roll_angle_right;
    joint_vel(8) = dleg_hip_yaw_angle_right;
    joint_vel(9) = dleg_hip_pitch_angle_right;
    joint_vel(10) = dleg_knee_angle_right;
    joint_vel(11) = 0.0;
    joint_vel(12) = 0.0;
    joint_vel(13) = 0.0;

    constexpr int float_q = 7;
    constexpr int float_qd = 7;
    q_.block<14, 1>(7, 0) = joint_angle; // q = [x, y, z, quaternion,joint]; quaternion = [x y z w]

    qd_.block<14, 1>(6, 0) = joint_vel; // joint:left right
    pinocchio::forwardKinematics(model_, data_, q_, qd_);
    pinocchio::framesForwardKinematics(model_, data_, q_);

    pinocchio::computeJointJacobians(model_, data_, q_); // joint jacobian
    pinocchio::getJointJacobian(model_, data_, rfoot_joint_id, pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED, J_r);
    pinocchio::getJointJacobian(model_, data_, lfoot_joint_id, pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED, J_l);

    // com calculate
    pinocchio::centerOfMass(model_, data_, q_, qd_, true);

    com_vel = data_.vcom[0];
    com_pos = data_.com[0];

    Lfoot_vel = pinocchio::getFrameVelocity(model_, data_, Lfoot_id, pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED).linear();
    Rfoot_vel = pinocchio::getFrameVelocity(model_, data_, Rfoot_id, pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED).linear();

    act.linkJacobian_right = J_r.block<3, 6>(0, 13);
    act.linkJacobian_left = J_l.block<3, 6>(0, 6);

    torso_pos = data_.oMf[torso_id].translation();
    Lfoot_pos = data_.oMf[Lfoot_id].translation();
    Rfoot_pos = data_.oMf[Rfoot_id].translation();

    /*foot pos wrt their hip*/ // for slip calcu usage
    /*calculate FK*/

    // compute com position, w.r.t. torso and ground calculated by qref
    Eigen::Vector3d com_wrt_rfoot;
    Eigen::Vector3d com_wrt_rground;
    com_wrt_rfoot = (com_pos - (data_.oMf[Rfoot_id].translation()));
    Eigen::Vector3d com_wrt_lfoot;
    Eigen::Vector3d com_wrt_lground;
    com_wrt_lfoot = (com_pos - (data_.oMf[Lfoot_id].translation()));

    // compute com velocity, w.r.t. torso and ground calculated by qref
    Eigen::Vector3d vcom_wrt_rfoot;
    Eigen::Vector3d vcom_wrt_rground;
    vcom_wrt_rfoot = (com_vel - (Rfoot_vel));
    Eigen::Vector3d vcom_wrt_lfoot;
    Eigen::Vector3d vcom_wrt_lground;
    vcom_wrt_lfoot = (com_vel - (Lfoot_vel));

    // transfer to ground frame
    Eigen::Matrix<double, 3, 3> torso_rotationmatrix;
    torso_rotationmatrix = RPY2ROT(torsoAngle[0], torsoAngle[1], 0.0);
    // torso_rotationmatrix.transposeInPlace();//trans to world frame
    com_wrt_rground = torso_rotationmatrix * com_wrt_rfoot;   // pos
    vcom_wrt_rground = torso_rotationmatrix * vcom_wrt_rfoot; // vel

    com_wrt_lground = torso_rotationmatrix * com_wrt_lfoot;   // pos
    vcom_wrt_lground = torso_rotationmatrix * vcom_wrt_lfoot; // vel

    /*foot pos wrt their hip*/ // for slip calcu usage
    /*calculate FK*/

    // slip state w.r.t. world frame
    act.slipfootPos_left = -com_wrt_lground;
    act.slipfootPos_right = -com_wrt_rground;
    act.slipfootVel_left = -vcom_wrt_lground;
    act.slipfootVel_right = -vcom_wrt_rground;

    double leglength_projecton_2_sagittal = sqrt(pow(act.slipfootPos_right[0], 2) + pow(act.slipfootPos_right[2], 2));
    act.slip_leg_length = {act.slipfootPos_right.norm(), act.slipfootPos_left.norm()};
    act.slip_leg_pitch_angle = {-atan(act.slipfootPos_right[0] / act.slipfootPos_right[2]), -atan(act.slipfootPos_left[0] / act.slipfootPos_left[2])};
    act.slip_leg_roll_angle = {atan(act.slipfootPos_right[1] / leglength_projecton_2_sagittal), atan(act.slipfootPos_left[1] / leglength_projecton_2_sagittal)};

    // std::cout << "legpitch" << act.slip_leg_pitch_angle << "\n"
    //           << std::endl;

    // std::cout << "legroll" << act.slip_leg_roll_angle << "\n"
    //           << std::endl;

    // slip state w.r.t. world frame
    act.slipState << act.slip_leg_length[0], act.slip_leg_roll_angle[0], act.slip_leg_pitch_angle[0], act.slip_leg_length[1], act.slip_leg_roll_angle[1], act.slip_leg_pitch_angle[1]; // right 3 left 3

    act.slipStateV << vcom_wrt_rground.dot(com_wrt_rground.normalized()), 0, 0, vcom_wrt_lground.dot(com_wrt_lground.normalized()), 0, 0; // right 3 left 3

    slipState_output = act.slipState;   // 3 right 3 left
    slipStateV_output = act.slipStateV; // 3 right 3 left
}

void fn_slip_trans::slip2robot(robot_para &ref, Eigen::Matrix<double, 3, 1> &torsoAngle)
{
    /*slip2robot*/ // IK problem
    // variable define
    double ref_right_slip_leglength = ref.slipState(0);
    double ref_right_slip_legrollangle = ref.slipState(1);
    double ref_right_slip_legpitchangle = ref.slipState(2);
    double dref_right_slip_leglength = ref.slipStateV(0);
    double dref_right_slip_legrollangle = ref.slipStateV(1);
    double dref_right_slip_legpitchangle = ref.slipStateV(2);
    double ref_left_slip_leglength = ref.slipState(3);
    double ref_left_slip_legrollangle = ref.slipState(4);
    double ref_left_slip_legpitchangle = ref.slipState(5);
    double dref_left_slip_leglength = ref.slipStateV(3);
    double dref_left_slip_legrollangle = ref.slipStateV(4);
    double dref_left_slip_legpitchangle = ref.slipStateV(5);

    double rz = ref_right_slip_leglength * cos(ref_right_slip_legrollangle) * cos(ref_right_slip_legpitchangle);
    double rx = ref_right_slip_leglength * cos(ref_right_slip_legrollangle) * sin(ref_right_slip_legpitchangle);
    double ry = ref_right_slip_leglength * sin(ref_right_slip_legrollangle);

    double lz = ref_left_slip_leglength * cos(ref_left_slip_legrollangle) * cos(ref_left_slip_legpitchangle);
    double lx = ref_left_slip_leglength * cos(ref_left_slip_legrollangle) * sin(ref_left_slip_legpitchangle);
    double ly = ref_left_slip_leglength * sin(ref_left_slip_legrollangle);

    Eigen::Matrix<double, 6, 1> footvel;

    footvel = slipVel_2_footVel(ref_right_slip_leglength, ref_right_slip_legrollangle,
                                ref_right_slip_legpitchangle, ref_left_slip_leglength,
                                ref_left_slip_legrollangle, ref_left_slip_legpitchangle,
                                dref_right_slip_leglength, dref_right_slip_legrollangle,
                                dref_right_slip_legpitchangle, dref_left_slip_leglength,
                                dref_left_slip_legrollangle, dref_left_slip_legpitchangle);

    // targetVelocity << val1, val2, val3, val4, val5, val6;

    ref.foot_pos_slip_right << rx, ry, -rz;
    ref.foot_pos_slip_left << lx, ly, -lz;

    ref.foot_vel_slip_right << footvel[0], footvel[1], footvel[2];
    ref.foot_vel_slip_left << footvel[3], footvel[4], footvel[5];

    /*through ik qref to calculate update end effector position, for further jacobian*/
    joint_angleIK << ref.jointAngle[6], ref.jointAngle[7], ref.jointAngle[8], ref.jointAngle[9], ref.jointAngle[10], ref.jointAngle[11], 0, ref.jointAngle[0], ref.jointAngle[1], ref.jointAngle[2], ref.jointAngle[3], ref.jointAngle[4], ref.jointAngle[5], 0;           // first left, second right
    joint_velIK << ref.jointAngleV[6], ref.jointAngleV[7], ref.jointAngleV[8], ref.jointAngleV[9], ref.jointAngleV[10], ref.jointAngleV[11], 0, ref.jointAngleV[0], ref.jointAngleV[1], ref.jointAngleV[2], ref.jointAngleV[3], ref.jointAngleV[4], ref.jointAngleV[5], 0; // first left, second right

    // qref set the kinematic valuables
    q_IK.block<14, 1>(7, 0) = joint_angleIK;
    qd_IK.block<14, 1>(6, 0) = joint_velIK;

    // calculate
    pinocchio::framesForwardKinematics(model_IK, data_IK, q_IK);

    // get
    ref.contactPointPos_left.block(0, 0, 3, 1) = data_IK.oMf[Lfoot_id].translation();
    ref.contactPointPos_right.block(0, 0, 3, 1) = data_IK.oMf[Rfoot_id].translation();

    // qref based jacobian
    pinocchio::computeJointJacobians(model_IK, data_IK, q_IK);

    // right ik jacobian
    pinocchio::Data::Matrix6x Jr_IK;
    Jr_IK.resize(6, 20);
    Jr_IK.setZero();
    // left ik jacobian
    pinocchio::Data::Matrix6x Jl_IK;
    Jl_IK.resize(6, 20);
    Jl_IK.setZero();

    // traditional jacobian 3 linear 3 rotate, not the spatial vector
    pinocchio::getJointJacobian(model_IK, data_IK, rfoot_joint_id, pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED, Jr_IK);
    pinocchio::getJointJacobian(model_IK, data_IK, lfoot_joint_id, pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED, Jl_IK);

    // set ik ref jacobian, this jacobian is not the real pos jacobian, rather a ik based qred jacobian
    ref.linkJacobian_right = Jr_IK.block<3, 6>(0, 13);
    ref.linkJacobian_left = Jl_IK.block<3, 6>(0, 6);

    /*IK mathamatics*/ // target set
    Eigen::VectorXd targetPosition = Eigen::VectorXd::Zero(6);
    Eigen::VectorXd targetVelocity = Eigen::VectorXd::Zero(6);
    Eigen::VectorXd targetyawVelocity = Eigen::VectorXd::Zero(6);

    // com pos w.r.t. torso, calculated by qref
    pinocchio::centerOfMass(model_IK, data_IK, q_IK, qd_IK, true);
    Eigen::Vector3d com_base_qref;
    Eigen::Vector3d comv_base_qref;
    com_base_qref = data_IK.com[0];
    comv_base_qref = data_IK.vcom[0];

    Eigen::Matrix<double, 3, 3> torso_rotationmatrix;
    torso_rotationmatrix = RPY2ROT(torsoAngle[0], torsoAngle[1], 0.0);

    targetPosition << torso_rotationmatrix.transpose() * (ref.foot_pos_slip_right) + com_base_qref, torso_rotationmatrix.transpose() * (ref.foot_pos_slip_left) + com_base_qref;
    targetVelocity << torso_rotationmatrix.transpose() * ref.foot_vel_slip_right, torso_rotationmatrix.transpose() * ref.foot_vel_slip_left;

    targetyawVelocity << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    // calculate ik result is qref
    Eigen::VectorXd angleAnkle = Eigen::VectorXd::Zero(2);
    ikcalculate.slipgetikResult(ref, targetPosition, targetVelocity, angleAnkle, targetyawVelocity, Time_step);

    ref.jointAngle[10] = torsoAngle[1] - ref.jointAngle[8] - ref.jointAngle[9];
    ref.jointAngle[4] = torsoAngle[1] - ref.jointAngle[2] - ref.jointAngle[3]; // ankle need parallel with ground
    ref.jointAngle[11] = -torsoAngle[0] - ref.jointAngle[6];
    ref.jointAngle[5] = -torsoAngle[0] - ref.jointAngle[0]; // ankle need parallel with ground

    jointAngel_output = ref.jointAngle;   // 6right 6left
    jointAngelV_output = ref.jointAngleV; // 6right 6left
}

void fn_slip_trans::calculate_slipstate()
{
    this->in_para.jointAngle = this->jointAngel_input;
    this->in_para.jointAngleV = this->jointAngelV_input;

    robot2slip(this->in_para, this->torso_input);
}

void fn_slip_trans::calculate_jointangle()
{
    this->in_para.slipState = this->slipState_input;
    this->in_para.slipStateV = this->slipStateV_input;

    slip2robot(this->in_para, this->torso_input);
}

Eigen::Matrix3d fn_slip_trans::RPY2ROT(double rollIn, double pitchIn, double yawIn)
{
    static Eigen::Matrix3d m;

    double sr = sin(rollIn);
    double cr = cos(rollIn);
    double sp = sin(pitchIn);
    double cp = cos(pitchIn);
    double sy = sin(yawIn);
    double cy = cos(yawIn);

    m(0, 0) = cy * cp;
    m(1, 0) = sy * cp;
    m(2, 0) = -sp;

    m(0, 1) = -sy * cr + cy * sp * sr;
    m(1, 1) = cy * cr + sy * sp * sr;
    m(2, 1) = cp * sr;

    m(0, 2) = sy * sr + cy * sp * cr;
    m(1, 2) = -cy * sr + sy * sp * cr;
    m(2, 2) = cp * cr;

    return m;
}

Eigen::Matrix<double, 6, 1> fn_slip_trans::slipVel_2_footVel(
    double &rleglength, double &rlegrollangle,
    double &rlegpitchangle, double &lleglength,
    double &llegrollangle, double &llegpitchangle,
    double &drleglength, double &drlegrollangle,
    double &drlegpitchangle, double &dlleglength,
    double &dllegrollangle, double &dllegpitchangle)
{
    double b_rfootvel_tmp;
    double c_rfootvel_tmp;
    double d_rfootvel_tmp;
    double e_rfootvel_tmp;
    double rfootvel_tmp;
    rfootvel_tmp = std::sin(rlegpitchangle);
    b_rfootvel_tmp = std::cos(rlegrollangle);
    c_rfootvel_tmp = std::sin(rlegrollangle);
    d_rfootvel_tmp = std::cos(rlegpitchangle);
    e_rfootvel_tmp = drlegrollangle * rleglength;
    Eigen::Matrix<double, 6, 1> v_out;
    // right
    v_out[0] = (drleglength * b_rfootvel_tmp * rfootvel_tmp -
                e_rfootvel_tmp * c_rfootvel_tmp * rfootvel_tmp) +
               drlegpitchangle * rleglength *
                   b_rfootvel_tmp * d_rfootvel_tmp;
    v_out[1] = drleglength * c_rfootvel_tmp +
               e_rfootvel_tmp * b_rfootvel_tmp;
    v_out[2] = (e_rfootvel_tmp * d_rfootvel_tmp * c_rfootvel_tmp -
                drleglength *
                    std::cos(rlegrollangle) * d_rfootvel_tmp) +
               drlegpitchangle * rleglength *
                   std::cos(rlegrollangle) * rfootvel_tmp;
    rfootvel_tmp = std::sin(llegpitchangle);
    b_rfootvel_tmp = std::cos(llegrollangle);
    c_rfootvel_tmp = std::sin(llegrollangle);
    d_rfootvel_tmp = std::cos(llegpitchangle);
    e_rfootvel_tmp = dllegrollangle * lleglength;
    // left
    v_out[3] = (dlleglength * b_rfootvel_tmp * rfootvel_tmp -
                e_rfootvel_tmp * c_rfootvel_tmp * rfootvel_tmp) +
               dllegpitchangle * lleglength *
                   b_rfootvel_tmp * d_rfootvel_tmp;
    v_out[4] = dlleglength * c_rfootvel_tmp +
               e_rfootvel_tmp * b_rfootvel_tmp;
    v_out[5] = (e_rfootvel_tmp * d_rfootvel_tmp * c_rfootvel_tmp -
                dlleglength *
                    std::cos(llegrollangle) * d_rfootvel_tmp) +
               dllegpitchangle * lleglength *
                   std::cos(llegrollangle) * rfootvel_tmp;

    return v_out;
}