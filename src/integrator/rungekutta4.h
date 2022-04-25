#pragma once
#include <valarray>

#include "integrator.h"


class RungeKutta4 : public Integrator
{
public:
    RungeKutta4(class Box *, double = 0.01);
    double next_step();

private:
    double dt2;
    std::valarray<double> r_old, v_old;
    std::valarray<double> K1x, K1v, K2x, K2v, K3x, K3v, K4x, K4v;
};
