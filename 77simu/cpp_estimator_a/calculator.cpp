#include <cmath>
#include <Eigen/LU>
#include "local_library.h"

using namespace std;
using namespace Eigen;
//定数
const double mu_const = 3.986004418e14; //GM_E m^3/s^2
const double J2_const = 1.082636e-3; //無次元 重力J2項
const double Earth_Radius = 6378136.6; //m

const double diff_t = propagation_step_time; //1step当たりの秒数

Vector3d position_differential(const Vector3d& velocity);
Vector3d velocity_differential(const Vector3d& position, const Vector3d& velocity, const Vector3d& acceleration, const double& Cd);

Matrix10d update_M_matrix(Matrix10d M, const Vector3d position, const Vector3d velocity, const Vector3d acceleration, const double Cd);
Matrix10d calculate_A_matrix(const Vector3d position, const Vector3d velocity, const Vector3d acceleration, const double Cd);
Matrix10d calculate_Q_matrix();

void Runge_kutta(vector<Vector3d, aligned_allocator<Vector3d>>& position_estimate, vector<Vector3d, aligned_allocator<Vector3d>>& velocity_estimate, vector<Vector3d, aligned_allocator<Vector3d>>& acceleration_estimate, vector<double>& Cd_estimate, Matrix10d &M, vector<Vector10d, aligned_allocator<Vector10d>>& M_store_vector)
{
    int iteration_number = round(1.0/diff_t); //1秒間あたりのstep数。微妙なずれを防ぐ。

    int i;

	Vector3d position = position_estimate.back();
	Vector3d velocity = velocity_estimate.back();
	Vector3d accelration = acceleration_estimate.back();
	double Cd = Cd_estimate.back();

    for(i = 0;i < iteration_number;++i){
        Vector3d k0 = position_differential(velocity);
        Vector3d l0 = velocity_differential(position, velocity, accelration, Cd);

        Vector3d tmp_position = position + k0*diff_t/2.0;
        Vector3d tmp_velocity = velocity + l0*diff_t/2.0;
        Vector3d k1 = position_differential(tmp_velocity);
        Vector3d l1 = velocity_differential(tmp_position, tmp_velocity, accelration, Cd);

        tmp_position = position + k1*diff_t/2.0;
        tmp_velocity = velocity + l1*diff_t/2.0;
        Vector3d k2 = position_differential(tmp_velocity);
        Vector3d l2 = velocity_differential(tmp_position, tmp_velocity, accelration, Cd);

        tmp_position = position + k2*diff_t;
        tmp_velocity = velocity + l2*diff_t;
        Vector3d k3 = position_differential(tmp_velocity);
        Vector3d l3 = velocity_differential(tmp_position, tmp_velocity, accelration, Cd);

        position += diff_t*(k0 + 2.0*k1 + 2.0*k2 + k3)/6.0;
        velocity += diff_t*(l0 + 2.0*l1 + 2.0*l2 + l3)/6.0;
        position_estimate.push_back(position);
        velocity_estimate.push_back(velocity);
		acceleration_estimate.push_back(accelration);
		Cd_estimate.push_back(Cd);

        M = update_M_matrix(M, position, velocity, accelration, Cd);
		Vector10d traceM = Matrix_trace_vector(M);
		M_store_vector.push_back(traceM);
    }
}

void Kalman_Filter(vector<Vector3d, aligned_allocator<Vector3d>>& position_estimate, vector<Vector3d, aligned_allocator<Vector3d>>& velocity_estimate, vector<Vector3d, aligned_allocator<Vector3d>>& acceleration_estimate, vector<double>& Cd_estimate, Matrix10d &M, vector<Vector10d, aligned_allocator<Vector10d>>& M_store_vector, const Vector3d true_position, const Vector3d true_velocity, vector<Vector3d, aligned_allocator<Vector3d>>& position_observed, vector<Vector3d, aligned_allocator<Vector3d>>& velocity_observed, vector<Vector10d, aligned_allocator<Vector10d>>& estimate_error)
{
	int i;

	//事前推定値
	Vector3d estimate_position = position_estimate.back();
	Vector3d estimate_velocity = velocity_estimate.back();
	Vector3d estimate_acceleration = acceleration_estimate.back();
	double Cd = Cd_estimate.back();
	
	//観測値
	Vector3d observed_position = true_position;
	Vector3d observed_velocity = true_velocity;

	normal_distribution<> position_dist(0.0, position_error);
	normal_distribution<> velocity_dist(0.0, velocity_error);

	for(i = 0;i < 3;++i){
		observed_position(i) += position_dist(mt);
	}
	for(i = 0;i < 3;++i){
		observed_velocity(i) += velocity_dist(mt);
	}
	
	Matrix<double, 6, 10> H = MatrixXd::Zero(6, 10); //観測行列
	for(i = 0;i < 6;++i){
		H(i,i) = 1.0;
	}

	Matrix<double, 10, 6> K; //カルマンゲイン
	Matrix<double, 6, 6> R = MatrixXd::Zero(6, 6); //観測誤差共分散行列
	for(i = 0;i < 3;++i){
		R(i,i) = position_error*position_error;
	}
	for(i = 3;i < 6;++i){
		R(i,i) = velocity_error*velocity_error;
	}
	
	Matrix<double, 6, 6> tmp = R + H*M*H.transpose();
	
	K = M*H.transpose()*tmp.inverse();
	
	Vector6d z_vector; //観測値
	z_vector << observed_position(0), observed_position(1), observed_position(2), observed_velocity(0), observed_velocity(1), observed_velocity(2);

	Vector10d x_predict; //事前推定値
	x_predict << estimate_position(0), estimate_position(1), estimate_position(2), estimate_velocity(0), estimate_velocity(1), estimate_velocity(2), estimate_acceleration(0), estimate_acceleration(1), estimate_acceleration(2), Cd;

	Vector10d x_update = x_predict + K*(z_vector - H*x_predict); //事後推定値
	M = (MatrixXd::Identity(10, 10) - K*H)*M;

	Vector3d update_position;
	update_position << x_update(0), x_update(1), x_update(2);
	Vector3d update_velocity;
	update_velocity << x_update(3), x_update(4), x_update(5);
	Vector3d update_acceleration;
	update_acceleration << x_update(6), x_update(7), x_update(8);
	Cd = x_update(9);

	//更新
	position_estimate.back() = update_position;
	velocity_estimate.back() = update_velocity;
	acceleration_estimate.back() = update_acceleration;
	Cd_estimate.back() = Cd;
	Vector10d traceM = Matrix_trace_vector(M);
	M_store_vector.back() = traceM;

	//ログ回り
	position_observed.push_back(observed_position);
	velocity_observed.push_back(observed_velocity);
	Vector10d x_error = K*(z_vector - H*x_predict);
	estimate_error.push_back(x_error);
	
	return;
}

Vector3d position_differential(const Vector3d& velocity)
{
    return velocity;
}

Vector3d velocity_differential(const Vector3d& position, const Vector3d& velocity, const Vector3d& acceleration, const double& Cd)
{
    double r = position.norm();
    double v = velocity.norm();

    double z = position(2);

    double ac_norm = -mu_const/position.squaredNorm(); //2体の重力項
    double tmp_J2_coefficient = 3.0/2.0*mu_const*J2_const*pow(Earth_Radius, 2.0)/pow(r, 4.0); //J2項の係数

    Vector3d all_acceleration = position/r;

    all_acceleration(0) *= ac_norm - tmp_J2_coefficient*(1.0 - 5.0*pow(z/r, 2.0));
    all_acceleration(1) *= ac_norm - tmp_J2_coefficient*(1.0 - 5.0*pow(z/r, 2.0));
    all_acceleration(2) *= ac_norm - tmp_J2_coefficient*(3.0 - 5.0*pow(z/r, 2.0));

    all_acceleration -= Cd*v*velocity; //-Cd*V^2*(Vi/V) 大気抵抗

	all_acceleration += acceleration; //残りの摂動要素

    return all_acceleration;
}

Matrix10d update_M_matrix(Matrix10d M, const Vector3d position, const Vector3d velocity, const Vector3d acceleration,  const double Cd)
{
    Matrix10d A = calculate_A_matrix(position, velocity, acceleration, Cd); //Jacobi行列
	Matrix10d Phi = MatrixXd::Identity(10,10) + diff_t * A;

	Matrix10d Q = calculate_Q_matrix();

	M = Phi*M*Phi.transpose() + diff_t*diff_t*Q;

	return M;
}

Matrix10d calculate_A_matrix(const Vector3d position, const Vector3d velocity, const Vector3d acceleration, const double Cd)
{
    double r = position.norm();
    double v = velocity.norm();

    double x = position(0); double y = position(1); double z = position(2);
    double vx = velocity(0); double vy = velocity(1); double vz = velocity(2);

    double J2_coefficient = 3.0/2.0*mu_const*J2_const*Earth_Radius*Earth_Radius; 

    Matrix10d A = MatrixXd::Zero(10,10);
    A(0,3) = 1.0; A(1,4) = 1.0; A(2,5) = 1.0;

    A(3,0) = 3.0*mu_const*x*x/pow(r, 5.0) - mu_const/pow(r, 3.0) - J2_coefficient*(1.0/pow(r, 5.0) - 5.0*(x*x + z*z)/pow(r, 7.0) + 35.0*x*x*z*z/pow(r, 9.0));
    A(3,1) = 3.0*mu_const*x*y/pow(r, 5.0) - J2_coefficient*( -5.0*x*y/pow(r, 7.0) + 35.0*x*y*z*z/pow(r, 9.0));
    A(3,2) = 3.0*mu_const*x*z/pow(r, 5.0) - J2_coefficient*( -15.0*x*z/pow(r, 7.0) + 35.0*x*z*z*z/pow(r, 9.0));

    A(4,0) = 3.0*mu_const*x*y/pow(r, 5.0) - J2_coefficient*( -5.0*x*y/pow(r, 7.0) + 35.0*x*y*z*z/pow(r, 9.0));
    A(4,1) = 3.0*mu_const*y*y/pow(r, 5.0) - mu_const/pow(r, 3.0) - J2_coefficient*(1.0/pow(r, 5.0) - 5.0*(y*y + z*z)/pow(r, 7.0) + 35.0*y*y*z*z/pow(r, 9.0));
    A(4,2) = 3.0*mu_const*y*z/pow(r, 5.0) - J2_coefficient*( -15.0*y*z/pow(r, 7.0) + 35.0*y*z*z*z/pow(r, 9.0));

    A(5,0) = 3.0*mu_const*x*z/pow(r, 5.0) - J2_coefficient*( -15.0*x*z/pow(r, 7.0) + 35.0*x*z*z*z/pow(r, 9.0));
    A(5,1) = 3.0*mu_const*y*z/pow(r, 5.0) - J2_coefficient*( -15.0*y*z/pow(r, 7.0) + 35.0*y*z*z*z/pow(r, 9.0));
    A(5,2) = 3.0*mu_const*z*z/pow(r, 5.0) - mu_const/pow(r, 3.0) - J2_coefficient*(3.0/pow(r, 5.0) - 30.0*z*z/pow(r, 7.0) + 35.0*pow(z, 4.0)/pow(r, 9.0));

    A(3,3) = -Cd*(vx*vx/v + v);    A(3,4) = -Cd*vx*vy/v;    A(3,5) = -Cd*vx*vz/v;
    A(4,3) = -Cd*vx*vy/v;    A(4,4) = -Cd*(vy*vy/v + v);    A(4,5) = -Cd*vy*vz/v;
    A(5,3) = -Cd*vx*vz/v;    A(5,4) = -Cd*vy*vz/v;    A(5,5) = -Cd*(vz*vz/v + v);

	A(3,6) = 1.0;	A(4,7) = 1.0;	A(5,8) = 1.0;

    A(3,9) = -v*vx;    A(4,9) = -v*vy;    A(5,9) = -v*vz;

    return A;
}

Matrix10d calculate_Q_matrix()
{
	int i;

	Matrix10d Q = MatrixXd::Zero(10, 10);
	for(i = 0;i < 3;++i){
		Q(i,i) = position_noise*position_noise;
	}
	for(i = 3;i < 6;++i){
		Q(i,i) = velocity_noise*velocity_noise;
	}
	for(i = 6;i < 9;++i){
		Q(i,i) = acceleration_noise*acceleration_noise;
	}
	Q(9,9) = Cd_noise*Cd_noise;

	return Q;
}