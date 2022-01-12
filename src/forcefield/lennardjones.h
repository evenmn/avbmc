#pragma once

#include <string>
#include <vector>
#include <valarray>
#include <memory>

#include "forcefield.h"


class LennardJones : public ForceField
{
public:
    //LennardJones(class System *);
    LennardJones(class System *, std::string);
    void read_param_file(std::string);
    void allocate_memory();
    void free_memory();
    void sort_params();
    
    double comp_twobody_par(int, int, std::valarray<double>, std::valarray<double> &, bool);
    //double comp_energy_mol(std::vector<class Particle>, class Molecule*);
    double comp_energy_par(std::vector<class Particle>, int);
    double comp_energy_par(std::vector<class Particle>, int, std::valarray<double> &, bool);
    ~LennardJones();

private:
    int nline;

    // vectors to store raw data from param file
    std::vector<double> sigma_vec, epsilon_vec, rc_vec;

    // matrices to store sorted params
    double **sigma_mat = nullptr;
    double **epsilon_mat = nullptr;
    double **rc_sqrd_mat = nullptr;
    double **shift_mat = nullptr;

    std::vector<std::vector<int> > neigh_lists;
};
