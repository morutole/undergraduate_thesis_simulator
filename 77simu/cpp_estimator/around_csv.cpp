#include <fstream>
#include <sstream>
#include <iomanip>
#include "local_library.h"

using namespace std;

void read_csv(const string file_name, vector<string>& header, vector<vector<double>>& true_vec)
{
    int i;

    ifstream ifs(file_name);
    string line;
    bool head = true;

    while(getline(ifs, line)){
        stringstream ss(line);
        string columns;

        if(head){
            while(getline(ss, columns, ',')){
                header.push_back(columns);
            }
            int n = header.size();
            true_vec.resize(n);

            head = false;
            
            continue;
        }

        for(i = 0;i < header.size();++i){
            getline(ss, columns, ',');
            if(!columns.empty()){
                double value = stod(columns, nullptr);
                true_vec.at(i).push_back(value);
            }else{ //欠損用
                double value = -1e10;
                true_vec.at(i).push_back(value);
            }
        }
    }

    return;
}

void to_csv(const string file_name, const vector<Vector3d, aligned_allocator<Vector3d>>& position_estimate, const vector<Vector3d, aligned_allocator<Vector3d>>& velocity_estimate, const vector<double>& Cd_estimate)
{
    ofstream ofs(file_name);
    int i,j;
    int n = position_estimate.size();
    for(i = 0;i < n;i += log_period){
        for(j = 0;j < 3;++j) ofs << setprecision(12) << position_estimate.at(i)(j) << ',';
        for(j = 0;j < 3;++j) ofs << setprecision(12) << velocity_estimate.at(i)(j) << ',';
        ofs << setprecision(12) << Cd_estimate.at(i) << ',';
        ofs << endl;

        cout << i << endl;
    }

    return;
}