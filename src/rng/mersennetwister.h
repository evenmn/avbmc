#pragma once
#include <vector>
#include <random>

#include "rng.h"


class MersenneTwister : public RandomNumberGenerator
{
public:
    MersenneTwister();
    int next_int(int) override;
    double next_double() override;
    double next_gaussian(double, double) override;
    int choice(std::vector<double>) override;
    ~MersenneTwister() = default;

private:
    std::mt19937 generator;
};
