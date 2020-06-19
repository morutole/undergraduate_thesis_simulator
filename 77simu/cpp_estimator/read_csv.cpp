#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "local_library.h"

using namespace std;

void read_csv(string file_name, vector<string>& header, vector<vector<double>>& true_vec)
{
    int i,j;

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
            }else{ //回避用
                double value = -1e10;
                true_vec.at(i).push_back(value);
            }
        }
    }

    return;
}