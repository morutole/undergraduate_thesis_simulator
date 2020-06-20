#ifndef _local_library_h_
#define _local_library_h_

#include <iostream>
#include <vector>
#include <Eigen/Core>
#include <random>
using namespace std;
using namespace Eigen;

//乱数extern 実体はmainにある
extern random_device seed_gen;
extern mt19937 mt;

//関数 //read_csv
void read_csv(string file_name, vector<string>& header, vector<vector<double>>& true_vec);

//value_io
bool pick_true_value(const vector<string>& header, const vector<vector<double>>& true_value, vector<Vector3d, aligned_allocator<Vector3d>>& position_true, vector<Vector3d, aligned_allocator<Vector3d>>& velocity_true);
void initialize_estimate(const vector<Vector3d, aligned_allocator<Vector3d>>& position_true, const vector<Vector3d, aligned_allocator<Vector3d>>& velocity_true,  vector<Vector3d, aligned_allocator<Vector3d>>& position_estimate, vector<Vector3d, aligned_allocator<Vector3d>>& velocity_estimate, vector<double>& Cd_estimate);

//calculator
void Runge_kutta(vector<Vector3d, aligned_allocator<Vector3d>>& position_estimate, vector<Vector3d, aligned_allocator<Vector3d>>& velocity_estimate, vector<double>& Cd_estimate);

#endif