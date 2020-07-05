#include "local_library.h"

using namespace std;
using namespace Eigen;

//乱数の実体
random_device seed_gen;
mt19937 mt(42); //全ての答え(seed固定)

int main() 
{
    vector<string> header;
    vector<vector<double>> true_value;
    read_csv(header, true_value);

    vector<Vector3d, aligned_allocator<Vector3d>> position_true; //地球中心慣性座標系
    vector<Vector3d, aligned_allocator<Vector3d>> velocity_true; //地球中心慣性座標系

    if(!pick_true_value(header, true_value, position_true, velocity_true)) return 0;

    //推定値
    vector<Vector3d, aligned_allocator<Vector3d>> position_estimate; 
    vector<Vector3d, aligned_allocator<Vector3d>> velocity_estimate;
    vector<Vector3d, aligned_allocator<Vector3d>> acceleration_estimate;
    vector<double> Cd_estimate; //抵抗係数推定

    Matrix10d M; //誤差共分散行列
    vector<Vector10d, aligned_allocator<Vector10d>> M_store_vector; //推定値共分散行列対角成分のみ

    initialize_estimate(position_true, velocity_true, position_estimate, velocity_estimate, acceleration_estimate, Cd_estimate, M, M_store_vector);

    //観測された値(観測誤差あり)
    vector<Vector3d, aligned_allocator<Vector3d>> position_observed; 
    vector<Vector3d, aligned_allocator<Vector3d>> velocity_observed;
    vector<Vector10d, aligned_allocator<Vector10d>> estimate_error; //修正値
    
    int i;
    double percent = 0.0;

    for(i = 0;i < position_true.size()-1;++i){
        Runge_kutta(position_estimate, velocity_estimate, acceleration_estimate, Cd_estimate, M, M_store_vector);

        if((i+1)%GPS_period == 0){
            Vector3d positon = position_true.at(i+1);
            Vector3d velocity =  velocity_true.at(i+1);
            
            Kalman_Filter(position_estimate, velocity_estimate, acceleration_estimate, Cd_estimate, M, M_store_vector, positon, velocity, position_observed, velocity_observed, estimate_error);
        }

        if(progress_percentage("simulation", i, position_true.size()-1, percent)) percent += 1.0;
    }
    
    to_csv(position_estimate, velocity_estimate, acceleration_estimate, Cd_estimate, M_store_vector, estimate_error, position_true, velocity_true, position_observed, velocity_observed);

    return 0;
}