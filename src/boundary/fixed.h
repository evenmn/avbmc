#pragma once
#include <cassert>
#include "boundary.h"

class Fixed : public Boundary
{
public:
    Fixed(class Box* box_in, const vec length);
    bool correct_position(mat &pos);
    bool correct_velocity(mat &vel);   // needed for reflective boundaries
    bool correct_distance(mat &dist);  // needed for periodic boundaries
    double comp_volume();

private:
    vec length;
    double volume;
};
