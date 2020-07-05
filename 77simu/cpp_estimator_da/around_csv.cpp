#include <fstream>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <direct.h>
#include "local_library.h"

using namespace std;

string make_log_dir();
void make_initial_condition_txt(const string log_dir_path);

void read_csv(vector<string>& header, vector<vector<double>>& true_vec)
{
    int i;

    ifstream ifs(input_csv_name);
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

void to_csv(const vector<Vector3d, aligned_allocator<Vector3d>>& position_estimate, const vector<Vector3d, aligned_allocator<Vector3d>>& velocity_estimate, vector<Vector3d, aligned_allocator<Vector3d>>& acceleration_estimate, const vector<double>& Cd_estimate, const vector<Vector10d, aligned_allocator<Vector10d>>& M_store_vector, const vector<Vector10d, aligned_allocator<Vector10d>>& estimate_error, const vector<Vector3d, aligned_allocator<Vector3d>>& position_true, const vector<Vector3d, aligned_allocator<Vector3d>>& velocity_true, const vector<Vector3d, aligned_allocator<Vector3d>>& position_observed, const vector<Vector3d, aligned_allocator<Vector3d>>& velocity_observed)
{
    string log_dir_path = make_log_dir();

    int i,j;

    double percent = 0.0;
    ofstream ofs(log_dir_path + "/" + output_estimate_csv_name);
    int n = position_estimate.size();
    for(i = 0;i < n;i += log_period){
        for(j = 0;j < 3;++j) ofs << setprecision(12) << position_estimate.at(i)(j) << ',';
        for(j = 0;j < 3;++j) ofs << setprecision(12) << velocity_estimate.at(i)(j) << ',';
        for(j = 0;j < 3;++j) ofs << setprecision(12) << acceleration_estimate.at(i)(j) << ',';
        ofs << setprecision(12) << Cd_estimate.at(i) << ',';
        ofs << endl;

        if(progress_percentage("csv_estimate", i, n, percent)) percent += 1.0;
    }

    percent = 0.0;
    ofstream ofs2(log_dir_path + "/" + output_error_csv_name);
    n = estimate_error.size();
    for(i = 0;i < n;++i){
        for(j = 0;j < 10;++j) ofs2 << setprecision(12) << estimate_error.at(i)(j) << ',';
        ofs2 << endl;

        if(progress_percentage("csv_error", i, n, percent)) percent += 1.0;
    }

    percent = 0.0;
    ofstream ofs3(log_dir_path + "/" + output_true_csv_name);
    n = position_true.size();

    for(i = 0;i < n;++i){
        for(j = 0;j < 3;++j) ofs3 << setprecision(12) << position_true.at(i)(j) << ',';
        for(j = 0;j < 3;++j) ofs3 << setprecision(12) << velocity_true.at(i)(j) << ',';
        ofs3 << endl;

        if(progress_percentage("csv_true", i, n, percent)) percent += 1.0;
    }

    percent = 0.0;
    ofstream ofs4(log_dir_path + "/" + output_observed_csv_name);
    n = position_observed.size() + 1; //observeは1つ少ない

    for(i = 0;i < n;++i){
        if(i == 0){
            for(j = 0;j < 3;++j) ofs4 << setprecision(12) << position_estimate.at(i)(j) << ',';
            for(j = 0;j < 3;++j) ofs4 << setprecision(12) << velocity_estimate.at(i)(j) << ',';
        }else{
            for(j = 0;j < 3;++j) ofs4 << setprecision(12) << position_observed.at(i-1)(j) << ','; //observeは1つ少ない
            for(j = 0;j < 3;++j) ofs4 << setprecision(12) << velocity_observed.at(i-1)(j) << ',';
        }
        ofs4 << endl;

        if(progress_percentage("csv_observe", i, n, percent)) percent += 1.0;
    }

    percent = 0.0;
    ofstream ofs5(log_dir_path + "/" + output_M_csv_name);
    n = M_store_vector.size();

    for(i = 0;i < n;++i){
        for(j = 0;j < 10;++j) ofs5 << setprecision(12) << M_store_vector.at(i)(j) << ',';
        ofs5 << endl;

        if(progress_percentage("csv_M", i, n, percent)) percent += 1.0;
    }

    make_initial_condition_txt(log_dir_path);

    return;
}

string make_log_dir()
{
    string log_dir_path = "log";

    _mkdir(log_dir_path.c_str()); //失敗しても問題無し

    //現在時刻を取得
    auto timer = time(NULL);
    auto now = localtime(&timer);
    char start_time_c[64];
    strftime(start_time_c, 64, "%y%m%d_%H%M%S", now);
    string time_info = start_time_c;
    log_dir_path += "/" + time_info;

    _mkdir(log_dir_path.c_str());

    return log_dir_path;
}

void make_initial_condition_txt(const string log_dir_path)
{
    ofstream ofs(log_dir_path + "/initial_condition.txt");
    ofs << "propagation step time[s]: " << propagation_step_time << endl;
    ofs << "log step time[s]: " << propagation_step_time*log_period << endl;
    ofs << "GPS period: " << GPS_period << endl;
    ofs << "xyz 位置観測標準偏差[m]: " << position_error << endl;
    ofs << "Vxyz 速度観測標準偏差[m/s]: " << velocity_error << endl;
    ofs << "xyz 外乱[m]: " << position_noise << endl;
    ofs << "Vxyz 外乱[m/s]: " << velocity_noise << endl;
    ofs << "axyz 外乱[m/s^2]: " << acceleration_noise << endl;
    ofs << "Cd 外乱[/m]: " << Cd_noise << endl;
    ofs << "初期推定大気抵抗[N]: " << initial_estimate_airdrag_force << endl;

    return;
}