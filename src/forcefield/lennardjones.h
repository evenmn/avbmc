#pragma once

#include <string>
#include <vector>
#include <valarray>
#include <memory>

#include "forcefield.h"


class LennardJones : public ForceField
{
public:
    LennardJones(class Box *, const std::string &);
    void read_param_file(const std::string &);
    void allocate_memory();
    void free_memory();
    void sort_params();
    
    double comp_twobody_par(int, int, std::valarray<double>, std::valarray<double> &, bool);
    //double comp_energy_par(int, std::valarray<double> &, bool);
    double comp_energy_par_neigh0_eng0(unsigned int, std::valarray<double> &, bool) override;
    double comp_energy_par_neigh1_eng0(unsigned int, std::valarray<double> &, bool) override;
    double comp_energy_par_neigh1_eng1(unsigned int, std::valarray<double> &, bool) override;
    ~LennardJones();

private:
    int neigh_id;

    // vectors to store raw data from param file
    std::vector<double> sigma_vec, epsilon_vec, rc_vec;

    // matrices to store sorted params
    double **sigma_mat = nullptr;
    double **epsilon_mat = nullptr;
    double **rc_sqrd_mat = nullptr;
    double **shift_mat = nullptr;

    //std::vector<std::vector<int> > neigh_lists;
};
