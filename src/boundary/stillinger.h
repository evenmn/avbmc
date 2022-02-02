#pragma once
#include <vector>
#include <string>
#include <valarray>

#include "boundary.h"


class Stillinger : public Boundary
{
public:
    Stillinger(class Box*, double = 2.0);
    void set_crit(std::string, std::string, double);
    void update();
    void check(int, std::valarray<int> &, std::valarray<int> &);
    bool correct_position();
    //double comp_volume();
    ~Stillinger();

private:
    unsigned int ntype;
    double r_csq, v_c, **r_csq_mat;
    std::vector<std::vector<int> > neigh_lists;
    std::valarray<int> in_cluster;
    std::valarray<int> checked;

};
