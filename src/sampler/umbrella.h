#pragma once
#include <functional>
#include <valarray>

#include "sampler.h"


class Umbrella : public Sampler
{
public:
    Umbrella(class Box*, std::function<double (int)>, int = 100);
    double w(int);

private:
    std::valarray<double> tabulate(std::function<double (int)>, int);
    std::valarray<double> tabulated;
    int maxpar;
    std::function <double (int)> f;
};
