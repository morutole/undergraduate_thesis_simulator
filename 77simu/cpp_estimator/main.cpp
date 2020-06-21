#include "local_library.h"

using namespace std;
using namespace Eigen;

//乱数の実体
random_device seed_gen;
mt19937 mt(seed_gen());

int main() 
{
    string input_csv = "observe.csv";

    vector<string> header;
    vector<vector<double>> true_value;
    read_csv(input_csv, header, true_value);

    vector<Vector3d, aligned_allocator<Vector3d>> position_true; //地球中心座標系
    vector<Vector3d, aligned_allocator<Vector3d>> velocity_true; //地球中心座標系

    if(!pick_true_value(header, true_value, position_true, velocity_true)) return 0;

    vector<Vector3d, aligned_allocator<Vector3d>> position_estimate; 
    vector<Vector3d, aligned_allocator<Vector3d>> velocity_estimate;
    vector<double> Cd_estimate; //抵抗係数推定

    vector<Matrix7d, aligned_allocator<Matrix7d>> M_store_vector; //推定値共分散行列

    initialize_estimate(position_true, velocity_true, position_estimate, velocity_estimate, Cd_estimate, M_store_vector);

    int i;

    for(i = 0;i < position_true.size()-1;++i){
        cout << "simulator " << i << endl; 
        Runge_kutta(position_estimate, velocity_estimate, Cd_estimate, M_store_vector);

        if((i+1)%GPS_period == 0){
            Vector3d positon = position_true.at(i+1);
            Vector3d velocity =  velocity_true.at(i+1);

            Kalman_Filter(position_estimate, velocity_estimate, Cd_estimate, M_store_vector, positon, velocity);
        }
    }

    string output_csv = "estimate.csv";
    to_csv(output_csv, position_estimate, velocity_estimate, Cd_estimate);

    return 0;
}