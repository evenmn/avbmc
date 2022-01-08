#pragma once
#include <vector>
#include <valarray>
#include <memory>

#include "boundary.h"


class Stillinger : public Boundary
{
public:
    Stillinger(class Box*, double = 2.0);
    Stillinger(std::shared_ptr<Box>, double = 2.0);
    void update();
    void check(int, std::valarray<int> &, std::valarray<int> &);
    bool correct_position();
    bool correct_velocity();   // needed for reflective boundaries
    bool correct_distance();  // needed for periodic boundaries
    double comp_volume();

private:
    double r_csq, v_c;
    std::vector<std::vector<int> > neigh_lists;
    std::valarray<int> in_cluster;
    std::valarray<int> checked;
};
