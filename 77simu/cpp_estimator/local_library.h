#ifndef _local_library_h_
#define _local_library_h_

#include <iostream>
#include <vector>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

//定数
const double mu_const = 3.986004418e14; //GM_E m^3/s^2
const double J2_const = 1.082636e-3; //無次元
const double Earth_Radius = 6378136.6; //m

//関数
void read_csv(string file_name, vector<string>& header, vector<vector<double>>& true_vec);

bool pick_true_value(const vector<string>& header, const vector<vector<double>>& true_value, vector<Vector3d, aligned_allocator<Vector3d>>& position_true, vector<Vector3d, aligned_allocator<Vector3d>>& velocity_true);
void initialize_estimate(const vector<Vector3d, aligned_allocator<Vector3d>>& position_true, const vector<Vector3d, aligned_allocator<Vector3d>>& velocity_true,  vector<Vector3d, aligned_allocator<Vector3d>>& position_estimate, vector<Vector3d, aligned_allocator<Vector3d>>& velocity_estimate);

#endif