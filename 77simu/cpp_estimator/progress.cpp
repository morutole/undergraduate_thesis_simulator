#include "local_library.h"

using namespace std;

bool progress_percentage(const string name, const int i, const int n, const double percent)
{
    int border = (double)n*percent*0.01;
    if(i >= border){
        cout << name << " " << percent << '%' << " " << i << " " << n << endl;
        return true;
    }else return false;
}