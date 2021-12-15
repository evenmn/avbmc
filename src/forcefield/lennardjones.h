#pragma once

#include <string>
#include <vector>

#include "forcefield.h"


class LennardJones : public ForceField
{
public:
    LennardJones(class Box* box_in);
    LennardJones(class Box* box_in, std::string params);
    void read_param_file(std::string params);
    void sort_params();
    //double eval_acc_element(const mat positions, const int i, const int j, rowvec& acc, const bool comp_energy);
    //double eval_acc_par(const mat positions, const int i, rowvec& acc, const bool comp_energy);
    //double eval_acc(const mat positions, mat& accs, vec& potengs, const bool comp_energy);
    //double comp_force_par(const rowvec pos, rowvec& acc);
    
    double comp_energy_mol(std::vector<class Particle *>, class Molecule*);
    double comp_energy_par(std::vector<class Particle *>, int);
    //double update_force_par(const mat positions, const int i);
    //double update_force_all();
    ~LennardJones();

    //std::vector<int> build_neigh_list(const mat positions, const int i);
    //void build_neigh_lists();

private:
    int nline;

    // vectors to store raw data from param file
    //std::vector<std::string> label1_vec, label2_vec;
    std::vector<double> sigma_vec, epsilon_vec, rc_vec;

    // matrices to store sorted params
    double **sigma_mat, **epsilon_mat, **rc_sqrd_mat, **shift_mat;

    std::vector<std::vector<int> > neigh_lists;
};
