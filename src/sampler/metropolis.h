#pragma once
#include "sampler.h"


class Metropolis : public Sampler
{
public:
    Metropolis(class Box*);
    double w(int);
};
