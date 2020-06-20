#ifndef _local_library_h_
#define _local_library_h_

#include <iostream>
#include <vector>
#include <Eigen/Core>
#include <random>
using namespace std;
using namespace Eigen;

using Matrix7d = Matrix<double, 7, 7>;

const double propagation_step_time = 1e-3; //1step当たりの時間をいくつにするか(ルンゲクッタ)
const int log_period = 10; //何stepごとにログをとるか propagation_step_time * log_periodで何sごとにログが出るか決まる。

const double position_error = 0.5; //位置に乗る観測誤差[m]
const double velocity_error = 5e-3; //速度に乗る観測誤差[m/s];

//乱数extern 実体はmainにある
extern random_device seed_gen;
extern mt19937 mt;

//関数 //around_csv
void read_csv(const string file_name, vector<string>& header, vector<vector<double>>& true_vec);
void to_csv(const string file_name, const vector<Vector3d, aligned_allocator<Vector3d>>& position_estimate, const vector<Vector3d, aligned_allocator<Vector3d>>& velocity_estimate);

//value_io
bool pick_true_value(const vector<string>& header, const vector<vector<double>>& true_value, vector<Vector3d, aligned_allocator<Vector3d>>& position_true, vector<Vector3d, aligned_allocator<Vector3d>>& velocity_true);
void initialize_estimate(const vector<Vector3d, aligned_allocator<Vector3d>>& position_true, const vector<Vector3d, aligned_allocator<Vector3d>>& velocity_true,  vector<Vector3d, aligned_allocator<Vector3d>>& position_estimate, vector<Vector3d, aligned_allocator<Vector3d>>& velocity_estimate, vector<double>& Cd_estimate, Matrix7d& M, vector<Matrix7d, aligned_allocator<Matrix7d>>& M_store_vector);

//calculator
void Runge_kutta(vector<Vector3d, aligned_allocator<Vector3d>>& position_estimate, vector<Vector3d, aligned_allocator<Vector3d>>& velocity_estimate, vector<double>& Cd_estimate);

#endif