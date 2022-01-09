#pragma once
#include <valarray>

#include "boundary.h"


class Fixed : public Boundary
{
public:
    Fixed(class Box *, std::valarray<double>);
    bool correct_position();
    //bool correct_velocity();   // needed for reflective boundaries
    //bool correct_distance();  // needed for periodic boundaries
    double comp_volume();

private:
    std::valarray<double> length;
    unsigned int boxdim;
    double volume;
};
