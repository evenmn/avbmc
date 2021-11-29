#pragma once
#include <cmath>
#include "sampler.h"

class Metropolis : public Sampler
{
public:
    Metropolis(class Box* box_in);
    double w(int npar);
};
