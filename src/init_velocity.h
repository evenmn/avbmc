#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <valarray>


class Velocity
{
public:
    Velocity() {};
    virtual std::vector<std::valarray<double> > get_velocity(unsigned int, unsigned int) = 0;
};

class Zero : public Velocity
{
public:
    Zero();
    std::vector<std::valarray<double> > get_velocity(unsigned int, unsigned int = 3) override;
};

class Gauss : public Velocity
{
public:
    Gauss(class RandomNumberGenerator* rng_in, double mean_in, double var_in);
    std::vector<std::valarray<double> > get_velocity(unsigned int, unsigned int = 3) override;
private:
    double mean, var;
    class RandomNumberGenerator* rng = nullptr;
};

class Temp : public Gauss
{
public:
    Temp(class RandomNumberGenerator* rng_in, double temp_in);
};
