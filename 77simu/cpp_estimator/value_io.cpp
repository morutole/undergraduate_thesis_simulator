#include "local_library.h"

using namespace std;
using namespace Eigen;

bool pick_true_value(const vector<string>& header, const vector<vector<double>>& true_value, vector<Vector3d, aligned_allocator<Vector3d>>& position_true, vector<Vector3d, aligned_allocator<Vector3d>>& velocity_true)
{
    int i,j;

    string position_word = "sat_position_i";
    string velocity_word = "sat_velocity_i";

    vector<int> position_index;
    vector<int> velocity_index;

    for(i = 0;i < header.size();++i){
        string str = header.at(i);

        if(str.substr(0, position_word.size()) == position_word) position_index.push_back(i);
        if(str.substr(0, velocity_word.size()) == velocity_word) velocity_index.push_back(i);
    }

    if(position_index.size() != 3 || velocity_index.size() != 3){
        cout << "index size is not 3." << endl;
        cout << "something is wrong with csv." << endl;
        return false;
    }

    for(i = 0;i < true_value.at(0).size();++i){
        Vector3d tmp;
        for(j = 0;j < 3;++j){
            tmp(j) = true_value.at(position_index.at(j)).at(i);
        }
        position_true.push_back(tmp);
    }

    for(i = 0;i < true_value.at(0).size();++i){
        Vector3d tmp;
        for(j = 0;j < 3;++j){
            tmp(j) = true_value.at(velocity_index.at(j)).at(i);
        }
        velocity_true.push_back(tmp);
    }

    return true;
}

void initialize_estimate(const vector<Vector3d, aligned_allocator<Vector3d>>& position_true, const vector<Vector3d, aligned_allocator<Vector3d>>& velocity_true,  vector<Vector3d, aligned_allocator<Vector3d>>& position_estimate, vector<Vector3d, aligned_allocator<Vector3d>>& velocity_estimate, vector<double>& Cd_estimate, Matrix7d& M, vector<Matrix7d, aligned_allocator<Matrix7d>>& M_store_vector)
{
    int i;
    //初期値決め　わざと誤差を結構入れておく。
    uniform_real_distribution<> position_dist(-2.0*position_error, 2.0*position_error);
    uniform_real_distribution<> velocity_dist(-2.0*velocity_error, 2.0*velocity_error);

    Vector3d position = position_true.front();
    for(i = 0;i < 3;++i){
        position(i) += position_dist(mt); 
    }
    position_estimate.push_back(position);

    Vector3d velocity = velocity_true.front();   
    for(i = 0;i < 3;++i){
        velocity(i) += velocity_dist(mt);
    }
    velocity_estimate.push_back(velocity);

    double initial_airdragforce = 1.3e-5;
    double initial_Cd = 1.3e-5/velocity.squaredNorm();
    Cd_estimate.push_back(initial_Cd);

    for(i = 0;i < 3;++i){
        M(i) = (10.0*position_error)*(10.0*position_error); //位置の誤差
    }
    for(i = 4;i < 6;++i){
        M(i) = (10.0*velocity_error)*(10.0*velocity_error); //速度の誤差
    }
    M(6) = (initial_Cd/10.0)*(initial_Cd/10.0);

    M_store_vector.push_back(M);

    return;
}