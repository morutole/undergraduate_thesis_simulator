#include "local_library.h"

using namespace std;
using namespace Eigen;

//乱数の実体
random_device seed_gen;
mt19937 mt(seed_gen());

int main() 
{
    string csv_file_name = "observe.csv";

    vector<string> header;
    vector<vector<double>> true_value;
    read_csv(csv_file_name, header, true_value);

    vector<Vector3d, aligned_allocator<Vector3d>> position_true; //地球中心座標系
    vector<Vector3d, aligned_allocator<Vector3d>> velocity_true; //地球中心座標系

    if(!pick_true_value(header, true_value, position_true, velocity_true)) return 0;

    vector<Vector3d, aligned_allocator<Vector3d>> position_estimate; 
    vector<Vector3d, aligned_allocator<Vector3d>> velocity_estimate;
    vector<double> Cd_estimate; //抵抗係数推定

    initialize_estimate(position_true, velocity_true, position_estimate, velocity_estimate, Cd_estimate);

    int i;

    for(i = 0;i < position_true.size()-1;++i){
        Runge_kutta(position_estimate, velocity_estimate, Cd_estimate);

    }

    return 0;
}
