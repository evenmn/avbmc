#pragma once
#include <functional>
#include <valarray>

#include "sampler.h"


class Umbrella : public Sampler
{
public:
    Umbrella(class System*, std::function<double (int)>, int = 100);
    Umbrella(class System*, std::valarray<double>);
    double w(int) override;

private:
    std::valarray<double> tabulate(std::function<double (int)>, int);
    std::valarray<double> tabulated;
    int maxpar;
    std::function <double (int)> f;
};
