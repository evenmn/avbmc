#pragma once
#include "forcefield.h"

class LennardJones : public ForceField
{
public:
    LennardJones(class Box* box_in);
    LennardJones(class Box* box_in, string params);
    void read_param_file(string params);
    void sort_params();
    double eval_acc_element(const mat positions, const int i, const int j, rowvec& acc, const bool comp_energy);
    double eval_acc_par(const mat positions, const int i, rowvec& acc, const bool comp_energy);
    double eval_acc(const mat positions, mat& accs, vec& potengs, const bool comp_energy);
    double comp_force_par(const rowvec pos, rowvec& acc);

    vector<int> build_neigh_list(const mat positions, const int i);
    void build_neigh_lists(const mat positions);

private:
    int nline;

    // vectors to store raw data from param file
    vector<string> chem_symbol1_vec;
    vector<string> chem_symbol2_vec;
    vector<double> sigma_vec;
    vector<double> epsilon_vec;
    vector<double> rc_vec;

    // matrices to store sorted params
    mat sigma_mat, epsilon_mat, rc_sqrd_mat;

    vector<vector<int>> neigh_lists;
    string label = "Lennard Jones";
};
