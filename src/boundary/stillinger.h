#pragma once
#include "boundary.h"

class Stillinger : public Boundary
{
public:
    Stillinger(class Box* box_in, double r_c_in=2.);
    bool correct_position(mat &pos);
    bool correct_velocity(mat &vel);   // needed for reflective boundaries
    bool correct_distance(mat &dist);  // needed for periodic boundaries
    double comp_volume();

private:
    double r_c, v_c;
};
