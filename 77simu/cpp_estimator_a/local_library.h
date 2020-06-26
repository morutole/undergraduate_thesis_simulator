#ifndef _local_library_h_
#define _local_library_h_

#include <iostream>
#include <vector>
#include <Eigen/Core>
#include <random>
using namespace std;
using namespace Eigen;

using Matrix10d = Matrix<double, 10, 10>;
using Vector6d = Matrix<double, 6, 1>;
using Vector10d = Matrix<double, 10, 1>;

//ファイル名
const string input_csv_name = "observe.csv";
const string output_estimate_csv_name = "estimate.csv";
const string output_error_csv_name = "error.csv";
const string output_true_csv_name = "true.csv";
const string output_observed_csv_name = "observed.csv";

const double propagation_step_time = 1e-3; //1step当たりの時間をいくつにするか(ルンゲクッタ)
const int log_period = 10; //何stepごとにログをとるか propagation_step_time * log_periodで何sごとにログが出るか決まる。
const int GPS_period = 10; //何sごとにGPSから[r,v]に関する情報を渡すか。

const double position_error = 0.5; //位置に乗る観測誤差[m]
const double velocity_error = 5e-2; //速度に乗る観測誤差[m/s];

const double position_noise = 1.0; //位置に乗るノイズ[m]
const double velocity_noise = 0.1; //速度に乗るノイズ[m/s]
const double acceleration_noise = 1e-2; //加速度に乗るノイズ[m/s^2] TBD
const double Cd_noise = 1e-13; //[/m]

const double initial_estimate_airdrag_force = 1.3e-5; //[N]
extern double Cd_error; //[/m] 普通のCdだけでなく、rho*Cd*A/m全体としてのCd;

//乱数extern 実体はmainにある
extern random_device seed_gen;
extern mt19937 mt;

//関数 //around_csv
void read_csv(vector<string>& header, vector<vector<double>>& true_vec);
void to_csv(const vector<Vector3d, aligned_allocator<Vector3d>>& position_estimate, const vector<Vector3d, aligned_allocator<Vector3d>>& velocity_estimate, vector<Vector3d, aligned_allocator<Vector3d>>& acceleration_estimate, const vector<double>& Cd_estimate, const vector<Vector10d, aligned_allocator<Vector10d>>& estimate_error, const vector<Vector3d, aligned_allocator<Vector3d>>& position_true, const vector<Vector3d, aligned_allocator<Vector3d>>& velocity_true, const vector<Vector3d, aligned_allocator<Vector3d>>& position_observed, const vector<Vector3d, aligned_allocator<Vector3d>>& velocity_observed);

//value_io
bool pick_true_value(const vector<string>& header, const vector<vector<double>>& true_value, vector<Vector3d, aligned_allocator<Vector3d>>& position_true, vector<Vector3d, aligned_allocator<Vector3d>>& velocity_true);
void initialize_estimate(const vector<Vector3d, aligned_allocator<Vector3d>>& position_true, const vector<Vector3d, aligned_allocator<Vector3d>>& velocity_true,  vector<Vector3d, aligned_allocator<Vector3d>>& position_estimate, vector<Vector3d, aligned_allocator<Vector3d>>& velocity_estimate, vector<Vector3d, aligned_allocator<Vector3d>>& acceleration_estimate, vector<double>& Cd_estimate, vector<Matrix10d, aligned_allocator<Matrix10d>>& M_store_vector);

//calculator
void Runge_kutta(vector<Vector3d, aligned_allocator<Vector3d>>& position_estimate, vector<Vector3d, aligned_allocator<Vector3d>>& velocity_estimate, vector<Vector3d, aligned_allocator<Vector3d>>& acceleration_estimate, vector<double>& Cd_estimate, vector<Matrix10d, aligned_allocator<Matrix10d>>& M_store_vector);
void Kalman_Filter(vector<Vector3d, aligned_allocator<Vector3d>>& position_estimate, vector<Vector3d, aligned_allocator<Vector3d>>& velocity_estimate, vector<Vector3d, aligned_allocator<Vector3d>>& acceleration_estimate, vector<double>& Cd_estimate, vector<Matrix10d, aligned_allocator<Matrix10d>>& M_store_vector, const Vector3d true_position, const Vector3d true_velocity, vector<Vector3d, aligned_allocator<Vector3d>>& position_observed, vector<Vector3d, aligned_allocator<Vector3d>>& velocity_observed, vector<Vector10d, aligned_allocator<Vector10d>>& estimate_error);

//progress
bool progress_percentage(const string name, const int i, const int n, const double percent);

#endif