#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <armadillo>

using namespace std;
using namespace arma;

class Velocity
{
public:
    Velocity(class Box* box_in);
    virtual mat get_velocity(const int npar, const int ndim) = 0;

protected:
    class Box* box = nullptr;
};

class Gauss : public Velocity
{
public:
    Gauss(class Box* box_in, const double mean_in, const double var_in);
    mat get_velocity(const int npar, const int ndim);
private:
    double mean, var;
};

class Temp : public Velocity
{
public:
    Temp(class Box* box_in, const double temp_in);
    mat get_velocity(const int npar, const int ndim);
private:
    double temp;
};
