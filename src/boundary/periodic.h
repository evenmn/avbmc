#pragma once
#include <iostream>
#include <valarray>

#include "boundary.h"


class Periodic : public Boundary
{
public:
    Periodic(class Box *, std::valarray<double>);
    void correct_position(unsigned int);
    void correct_distance(std::valarray<double> &);   // needed for periodic boundaries

private:
    unsigned int ndim;
    double volume;
    std::valarray<double> length;
};
