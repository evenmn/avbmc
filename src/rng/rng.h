#pragma once
#include <iostream>
#include <random>
#include <vector>
#include <algorithm>
#include <cassert>


class RandomNumberGenerator
{
public:
    RandomNumberGenerator() : generator(seed()) {}
    virtual int next_int(int upper_limit) = 0;
    virtual double next_double() = 0;
    virtual double next_gaussian(double mean=0., double variance=1.) = 0;
    virtual int choice(std::vector<double> probabilities) = 0;
    virtual ~RandomNumberGenerator() = default;

    std::random_device seed;
    std::mt19937 generator;
};
