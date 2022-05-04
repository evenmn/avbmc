#pragma once
#include <random>
#include <vector>
#include <string>


class RandomNumberGenerator
{
public:
    RandomNumberGenerator() {};
    virtual void set_seed(unsigned int) = 0;
    virtual int next_int(int) = 0;
    virtual double next_double() = 0;
    virtual double next_gaussian(double = 0.0, double = 1.0) = 0;
    virtual int choice(std::vector<double>) = 0;
    virtual ~RandomNumberGenerator() = default;

    std::string label;

protected:
    std::random_device seed;
};
