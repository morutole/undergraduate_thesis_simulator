#include <iostream>
#include <Eigen/Dense>
#include "local_library.h"
using namespace std;
using namespace Eigen;

int main() 
{
    string csv_file_name = "observe.csv";

    vector<string> header;
    vector<vector<double>> true_value;
    read_csv(csv_file_name, header, true_value);

    cout << csv_file_name.size() << endl;
    cout << true_value.at(0).size() << endl;

    return 0;
}
