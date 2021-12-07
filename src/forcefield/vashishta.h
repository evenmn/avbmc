#pragma once
#include <stdlib.h>
#include "forcefield.h"

class Vashishta : public ForceField
{
public:
    Vashishta(class Box* box_in, std::string params);
    void read_param_file(std::string params);
    void sort_params();
    //double eval_acc_element(const mat positions, const int i, const int j, rowvec& acc, const bool comp_energy);
    //double eval_acc_par(const mat positions, const int i, rowvec& acc, const bool comp_energy);
    //double eval_acc(const mat positions, mat& accs, vec& potengs, const bool comp_energy);
    //double comp_force_par(const rowvec pos, rowvec& acc);
    
    double comp_twobody_par(std::vector<class Particle *>, int);
    double comp_threebody_par(std::vector<class Particle *>, int, int);

    double comp_energy_mol(std::vector<class Particle *>, class Molecule);
    double comp_energy_par(std::vector<class Particle *>, class Particle* atom);
    //double update_force_par(const mat positions, const int i);
    //double update_force_all();
    ~Vashishta();

    //std::vector<int> build_neigh_list(const mat positions, const int i);
    //void build_neigh_lists();

private:
    int nline;

    // vectors to store raw data from param file
    std::vector<std::string> label1_vec, label2_vec, label3_vec;
    std::vector<double> H_vec, eta_vec, Zi_vec, Zj_vec, lambda1_vec, D_vec, lambda4_vec;
    std::vector<double> W_vec, rc_vec, B_vec, gamma_vec, r0_vec, C_vec, costheta_vec;

    // matrices to store sorted params
    double ***H_mat, ***eta_mat, ***Zi_mat, ***Zj_mat, ***lambda1_mat, ***D_mat, ***lambda4_mat;
    double ***W_mat, ***rc_sqrd_mat, ***B_mat, ***gamma_mat, ***r0_mat, ***C_mat, ***costheta_mat;

    std::vector<std::vector<int> > neigh_lists;
};
