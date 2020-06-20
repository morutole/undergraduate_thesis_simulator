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

void initialize_estimate(const vector<Vector3d, aligned_allocator<Vector3d>>& position_true, const vector<Vector3d, aligned_allocator<Vector3d>>& velocity_true,  vector<Vector3d, aligned_allocator<Vector3d>>& position_estimate, vector<Vector3d, aligned_allocator<Vector3d>>& velocity_estimate, vector<double>& Cd_estimate)
{
    int i;
    //初期値決め
    uniform_real_distribution<> position_dist(-1.0, 1.0); //精度50cm
    uniform_real_distribution<> velocity_dist(-0.01, 0.01); //精度1cm/s

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

    double initial_airdragforce = 1.3e5;
    double initial_Cd = 1.3e5/velocity.squaredNorm();
    Cd_estimate.push_back(initial_Cd);

    return;
}