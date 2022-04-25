#pragma once
#include <vector>
#include <string>
#include <valarray>

#include "constraint.h"


class Stillinger : public Constraint
{
public:
    Stillinger(class Box*, double = 2.0);
    void set_criterion(std::string, std::string, double);
    void check_neigh_recu(int, std::valarray<char> &, std::valarray<char> &);
    bool verify();
    //double comp_volume();
    ~Stillinger();

private:
    unsigned int ntype, cutoff_id;
    double v_c, **r_csq_mat;
    std::vector<std::vector<int> > neigh_lists;
    std::valarray<char> in_cluster, checked;
};
