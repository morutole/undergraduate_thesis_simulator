#include <cmath>
#include "local_library.h"

using namespace std;
using namespace Eigen;
//定数
const double mu_const = 3.986004418e14; //GM_E m^3/s^2
const double J2_const = 1.082636e-3; //無次元
const double Earth_Radius = 6378136.6; //m

Vector3d position_differential(const Vector3d& velocity);
Vector3d velocity_differential(const Vector3d& position, const Vector3d& velocity, const double& Cd);

void Runge_kutta(vector<Vector3d, aligned_allocator<Vector3d>>& position_estimate, vector<Vector3d, aligned_allocator<Vector3d>>& velocity_estimate, vector<double>& Cd_estimate)
{
    double diff_t = 1e-3; //1step当たりの時間をいくつにするか
    int iteration_number = round(1.0/diff_t); //1秒間あたりのstep数。微妙なずれを防ぐ。

    int i;

    for(i = 0;i < iteration_number;++i){
        Vector3d position = position_estimate.back();
        Vector3d velocity = velocity_estimate.back();
        double Cd = Cd_estimate.back();

        Vector3d k0 = position_differential(velocity);
        Vector3d l0 = velocity_differential(position, velocity, Cd);

        Vector3d tmp_position = position + k0*diff_t/2.0;
        Vector3d tmp_velocity = velocity + l0*diff_t/2.0;
        Vector3d k1 = position_differential(tmp_velocity);
        Vector3d l1 = velocity_differential(tmp_position, tmp_velocity, Cd);

        tmp_position = position + k1*diff_t/2.0;
        tmp_velocity = velocity + l1*diff_t/2.0;
        Vector3d k2 = position_differential(tmp_velocity);
        Vector3d l2 = velocity_differential(tmp_position, tmp_velocity, Cd);

        tmp_position = position + k2*diff_t;
        tmp_velocity = velocity + l2*diff_t;
        Vector3d k3 = position_differential(tmp_velocity);
        Vector3d l3 = velocity_differential(tmp_position, tmp_velocity, Cd);

        position += diff_t*(k0 + 2.0*k1 + 2.0*k2 + k3)/6.0;
        velocity += diff_t*(l0 + 2.0*l1 + 2.0*l2 + l3)/6.0;
        position_estimate.push_back(position);
        velocity_estimate.push_back(velocity);
    }
}

Vector3d position_differential(const Vector3d& velocity)
{
    return velocity;
}

Vector3d velocity_differential(const Vector3d& position, const Vector3d& velocity, const double& Cd)
{
    double r = position.norm();
    double v = velocity.norm();

    double z = position(2);

    double ac_norm = -mu_const/position.squaredNorm(); //2体の重力項
    double tmp_J2_coefficient = 3.0/2.0*mu_const*J2_const*pow(Earth_Radius, 2.0)/pow(r, 4.0); //J2項の係数

    Vector3d acceleration = position/r;

    acceleration(0) *= ac_norm + tmp_J2_coefficient*(1.0 - 5.0*pow(z/r, 2.0));
    acceleration(1) *= ac_norm + tmp_J2_coefficient*(1.0 - 5.0*pow(z/r, 2.0));
    acceleration(2) *= ac_norm + tmp_J2_coefficient*(3.0 - 5.0*pow(z/r, 2.0));

    acceleration -= Cd*v*velocity; //-Cd*V^2*(Vi/V) 大気抵抗

    return acceleration;
}