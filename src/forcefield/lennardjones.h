#pragma once
#include "forcefield.h"

class LennardJones : public ForceField
{
public:
    LennardJones(class Box* box_in, string params);
    void read_param_file(string params);
    double eval_acc_par(const mat positions, const int i, rowvec& acc, const bool comp_energy=false);
    double eval_acc(const mat positions, mat& accs, const bool comp_energy=false);

    vector<int> build_neigh_list(const mat positions, const int i);
    void build_neigh_lists(const mat positions);

private:
    mat sigma_mat, epsilon_mat, rc_sqrd_mat;

    vector<vector<int>> neigh_lists;
    string label = "Lennard Jones";
};
