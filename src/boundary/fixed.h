#pragma once
#include <valarray>

#include "boundary.h"


class Fixed : public Boundary
{
public:
    Fixed(class Box *, std::valarray<double>);
    bool correct_position();
    double comp_volume();

private:
    std::valarray<double> length;
    unsigned int boxdim;
    double volume;
};
