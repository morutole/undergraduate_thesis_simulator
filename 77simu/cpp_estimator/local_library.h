#ifndef _local_library_h_
#define _local_library_h_

#include <iostream>
#include <vector>
using namespace std;

void read_csv(string file_name, vector<string>& header, vector<vector<double>>& true_vec);

#endif