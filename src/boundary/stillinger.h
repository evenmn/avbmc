#pragma once
#include <string>
#include <vector>
#include <valarray>

#include "boundary.h"

class Stillinger : public Boundary
{
public:
    Stillinger(class Box* box_in, double r_c_in=2.);
    void update();
    void check(const int i, std::valarray<bool> &in_cluster, std::valarray<bool> &checked);
    bool correct_position();
    bool correct_velocity();   // needed for reflective boundaries
    bool correct_distance();  // needed for periodic boundaries
    double comp_volume();

private:
    double r_c, v_c;

    std::vector<std::vector<int> > neigh_lists;
};
