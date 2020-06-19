#include <random>
#include "local_library.h"

using namespace std;
using namespace Eigen;

//乱数
random_device seed_gen;
mt19937 mt(seed_gen());

bool pick_true_value(const vector<string>& header, const vector<vector<double>>& true_value, vector<Vector3d, aligned_allocator<Vector3d>>& position_true, vector<Vector3d, aligned_allocator<Vector3d>>& velocity_true)
{
    int i;

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

    return true;
}

void initialize_estimate(const vector<Vector3d, aligned_allocator<Vector3d>>& position_true, const vector<Vector3d, aligned_allocator<Vector3d>>& velocity_true,  vector<Vector3d, aligned_allocator<Vector3d>>& position_estimate, vector<Vector3d, aligned_allocator<Vector3d>>& velocity_estimate)
{




    return;
}