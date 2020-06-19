#include <iostream>
#include <Eigen/Dense>
#include "local_library.h"
using namespace std;
using namespace Eigen;

int main() 
{
    cout << "Hello world" << endl;
    Vector3d tmp;
    tmp << 3, 4, 5;
    cout << tmp.norm() << endl;
    good_bye();
    return 0;
}
