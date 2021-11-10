#pragma once
#include <functional>
#include <valarray>

#include "sampler.h"

class Umbrella : public Sampler
{
public:
    Umbrella(class Box* box_in, std::function<double (int)> f_in, const int maxpar=100);

    double get_acceptance_prob(class Moves* move, const double temp, const double chempot);

private:
    std::valarray<double> tabulate(std::function<double (int)> f_in, const int maxpar);

    std::function<double (int)> f;
    std::valarray<double> tabulated;
    int maxpar;
};