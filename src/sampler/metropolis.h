#pragma once
#include "sampler.h"


class Metropolis : public Sampler
{
public:
    Metropolis(class System*);
    double w(int);
};
