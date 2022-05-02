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
    Velocity();
    virtual mat get_velocity(const int npar, const int ndim) = 0;
};

class Zero : public Velocity
{
public:
    Zero();
    mat get_velocity(const int npar, const int ndim) override;
};

class Gauss : public Velocity
{
public:
    Gauss(class RandomNumberGenerator* rng_in, const double mean_in, const double var_in);
    mat get_velocity(const int npar, const int ndim) override;
private:
    double mean, var;
    class RandomNumberGenerator* rng = nullptr;
};

class Temp : public Gauss
{
public:
    Temp(class RandomNumberGenerator* rng_in, const double temp_in);
};
