#pragma once
#include <vector>
#include <random>

#include "rng.h"


class MersenneTwister : public RandomNumberGenerator
{
public:
    MersenneTwister();
    int next_int(int);
    double next_double();
    double next_gaussian(double, double);
    int choice(std::vector<double>);

private:
    std::mt19937 generator;
};
