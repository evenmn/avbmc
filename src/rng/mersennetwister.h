#pragma once
#include "rng.h"

class MersenneTwister : public RandomNumberGenerator
{
public:
    MersenneTwister();
    int next_int(int upper_limit);
    double next_double();
    double next_gaussian(double mean, double variance);
    int choice(std::vector<double> probabilities);
};
