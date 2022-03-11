#pragma once
#include <string>
#include <vector>
#include <valarray>

#include "forcefield.h"


class Vashishta : public ForceField
{
public:
    Vashishta(class Box *, std::string);
    void read_param_file(std::string);
    void allocate_memory();
    void free_memory();
    void sort_params();
    
    double comp_twobody_par(int, int, double, std::valarray<double> &, bool);
    double comp_threebody_par(int, int, int, std::valarray<double>, std::valarray<double>, double, std::valarray<double> &, bool);
    double comp_threebody_par(int, int, int, double, double, std::valarray<double>, std::valarray<double>, std::valarray<double> &, bool);
    double comp_energy_par(int, std::valarray<double> &, bool);
    double comp_energy_par_neigh(int, std::valarray<double> &, bool);
    double comp_energy_par_noneigh(int, std::valarray<double> &, bool);

    ~Vashishta();

private:
    // vectors to store raw data from param file
    unsigned int neigh_id_rc, neigh_id_r0;
    std::vector<std::string> label3_vec;
    std::vector<double> H_vec, eta_vec, Zi_vec, Zj_vec, lambda1_vec, D_vec, lambda4_vec;
    std::vector<double> W_vec, rc_vec, B_vec, gamma_vec, r0_vec, C_vec, costheta_vec;

    // matrices to store sorted params
    double **H_mat, **eta_mat, **Zi_mat, **Zj_mat, **lambda1inv_mat, **D_mat, **lambda4inv_mat;
    double **W_mat, **rc_mat, ***B_mat, **gamma_mat, **r0_mat, **shift_mat, ***C_mat, ***costheta_mat;

    std::vector<std::vector<int> > neigh_lists;
};
