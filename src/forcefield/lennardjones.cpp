#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cassert>

#include "lennardjones.h"
#include "../box.h"
#include "../system.h"
#include "../particle.h"
#include "../distance_manager.h"


/* ----------------------------------------------------------------------------
   This is the default constructor when a parameter file 'params' is given.
------------------------------------------------------------------------------- */

LennardJones::LennardJones(System* system_in, const std::string params)
    : ForceField(system_in)
{
    label = "Lennard-Jones";
    paramfile = params;
    read_param_file(params);
    create_label_mapping();
    allocate_memory();
    sort_params();
}


/* ----------------------------------------------------------------------------
   Read parameter file 'params' and store parameters globally. It takes the
   following form:
       <label1> <label2> <sigma> <epsilon> <rc>
------------------------------------------------------------------------------- */

void LennardJones::read_param_file(const std::string params)
{
    nline = 0;
    std::ifstream infile(params);
    if (infile.is_open()) {
        std::string line, label1, label2;
        double sigma, epsilon, rc;
        while (std::getline(infile, line))
        {
            std::istringstream iss(line);
            if (line.rfind("#", 0) == 0) {
                // comments are allowed in parameter file
                continue;
            }
            else if (line.empty()) {
                // empty lines are allowed in parameter file
                continue;
            }
            else if (iss >> label1 >> label2 >> sigma >> epsilon >> rc) { 
                label1_vec.push_back(label1);
                label2_vec.push_back(label2);
                sigma_vec.push_back(sigma);
                epsilon_vec.push_back(epsilon);
                rc_vec.push_back(rc);
                nline ++;
            }
            else {
                std::cout << "Warning: Corrupt line in parameter file!" << std::endl;
                std::cout << "Ignoring line: '" + line + "'" << std::endl;
            }
        }
    }
    else {
        std::cout << "\nParameter file '" + params + "' was not found!" << std::endl;
        exit(0);
    }
}


/* ----------------------------------------------------------------------------
   Parameters have to be sorted with respect to the particle types, and are
   stored in matrices.
------------------------------------------------------------------------------- */

void LennardJones::sort_params()
{
    // link list of labels to list of type indices
    double rcsq, s6, s12;
    unsigned int type1, type2, i;
    std::vector<int> types1_vec, types2_vec;
    for(std::string label : label1_vec){
        types1_vec.push_back(label2type.at(label));
    }
    
    for(std::string label : label2_vec){
        types2_vec.push_back(label2type.at(label));
    }

    // fill up matrices with parameters
    for (i=0; i<nline; i++) {
        type1 = types1_vec[i];
        type2 = types2_vec[i];
        sigma_mat[type1][type2] = sigma_vec[i];
        sigma_mat[type2][type1] = sigma_vec[i];
        epsilon_mat[type1][type2] = epsilon_vec[i];
        epsilon_mat[type2][type1] = epsilon_vec[i];
        rcsq = rc_vec[i] * rc_vec[i];
        rc_sqrd_mat[type1][type2] = rcsq;
        rc_sqrd_mat[type2][type1] = rcsq;
        s6 = std::pow(sigma_vec[i] / rcsq, 3);
        s12 = s6 * s6;
        shift_mat[type1][type2] = s12 - s6;
        shift_mat[type2][type1] = s12 - s6;
    }
}


/* ----------------------------------------------------------------------------
   Set ids of neighbor lists associated with energy cutoff. This has to be
   done after the boxes are defined, but before the simulation is started.
------------------------------------------------------------------------------- */

void LennardJones::set_cutoff_ids()
{
    for (Box* box : system->boxes) {
        neigh_ids.push_back(box->distance_manager->add_cutoff(rc_sqrd_mat));
    }
}


/* ----------------------------------------------------------------------------
   Compute interaction energy between two particles of types 'typei' and
   'typej', respectively, separated by a distance vector 'delij'. Updates a
   force array 'force' if 'comp_force' is true.
------------------------------------------------------------------------------- */

double LennardJones::comp_twobody_par(const int typei, const int typej,
                                      const std::valarray<double> delij,
                                      std::valarray<double> &force, const bool comp_force)
{
    double rijsq, rijinvsq, s6, s12, energy;

    energy = 0.;
    rijsq = norm(delij); 
    if (rijsq < rc_sqrd_mat[typei][typej]) {
        rijinvsq = 1. / rijsq;
        s6 = std::pow(sigma_mat[typei][typej] * rijinvsq, 3);
        s12 = s6 * s6;
        energy = epsilon_mat[typei][typej] * (s12 - s6 - shift_mat[typei][typej]);
        if (comp_force) {
            force += epsilon_mat[typei][typej] * (2. * s6 - s12) * delij * rijinvsq;
        }
    }
    return energy;
}


/* ----------------------------------------------------------------------------
   Compute energy contribution from a particle 'i'
------------------------------------------------------------------------------- */

double LennardJones::comp_energy_par_noneigh(Box* box, const int i,
                                     std::valarray<double> &force, const bool comp_force)
{
    // declare variables
    int npar, typei, typej, j;
    npar = box->particles.size();
    typei = box->particles[i].type;
    std::valarray<double> delij;
    force.resize(system->ndim, 0.);

    double energy = 0.;
    for (j=0; j<i; j++) {
        typej = box->particles[j].type;
        delij = box->particles[j].r - box->particles[i].r;
        energy += comp_twobody_par(typei, typej, delij, force, comp_force);
    }
    for (j=i+1; j<npar; j++) {
        typej = box->particles[j].type;
        delij = box->particles[j].r - box->particles[i].r;
        energy += comp_twobody_par(typei, typej, delij, force, comp_force);
    }
    force *= 24;
    return 4 * energy;
}


/* ----------------------------------------------------------------------------
   Compute energy contribution from a particle 'i'
------------------------------------------------------------------------------- */

double LennardJones::comp_energy_par_neigh(Box* box, const int i, std::valarray<double> &force,
                                     const bool comp_force)
{
    // declare variables
    int npar, typei, typej;
    double rijsq, rijinvsq, s6, s12, energy, energyij;

    npar = box->particles.size();
    typei = box->particles[i].type;
    std::valarray<double> delij;
    force.resize(system->ndim, 0.);
    std::vector<std::vector<int> > neigh_list;

    neigh_list = box->distance_manager->neigh_lists[neigh_ids[box->box_id]];

    energy = 0.;
    for (int j : neigh_list[i]) {
        typej = box->particles[j].type;
        rijsq = box->distance_manager->distance_mat[i][j];
        rijinvsq = 1. / rijsq;
        s6 = std::pow(sigma_mat[typei][typej] * rijinvsq, 3);
        s12 = s6 * s6;
        energyij =  epsilon_mat[typei][typej] * (s12 - s6 - shift_mat[typei][typej]);
        if (comp_force) {
            delij = box->distance_manager->distance_cube[i][j];
            force += epsilon_mat[typei][typej] * (2. * s6 - s12) * delij * rijinvsq;
        }
        //if (system->store_energies) {
        //    poteng_mat[i][j] = energyij;
        //    poteng_mat[j][i] = energyij;
        //}
        energy += energyij;
    }
    //if (system->store_energies) {
    //    poteng[i] = energy;
    //}
    force *= 24;
    return 4 * energy;
}


/* ----------------------------------------------------------------------------
   Forwarding the computations of the energy of a particle 'i' to other 
   functions
------------------------------------------------------------------------------- */

double LennardJones::comp_energy_par(Box* box, const int i, std::valarray<double> &force,
                                const bool comp_force)
{
    double energy;

    if (system->store_distances) {
        comp_energy_par_neigh(box, i, force, comp_force);
    }
    else {
        comp_energy_par_noneigh(box, i, force, comp_force);
    }
    return energy;
}


/* ----------------------------------------------------------------------------
   Allocate memory for matrices
------------------------------------------------------------------------------- */

void LennardJones::allocate_memory()
{
    sigma_mat = new double*[ntype];
    epsilon_mat = new double*[ntype];
    rc_sqrd_mat = new double*[ntype];
    shift_mat = new double*[ntype];
    for (unsigned int i=0; i<ntype; i++) {
        sigma_mat[i] = new double[ntype];
        epsilon_mat[i] = new double[ntype];
        rc_sqrd_mat[i] = new double[ntype];
        shift_mat[i] = new double[ntype];
    }
}


/* ----------------------------------------------------------------------------
   Free memory for matrices
------------------------------------------------------------------------------- */

void LennardJones::free_memory()
{
    for (unsigned int i = 0; i < ntype; i++) {
        delete[] sigma_mat[i];
        delete[] epsilon_mat[i];
        delete[] rc_sqrd_mat[i];
        delete[] shift_mat[i];
    }
    delete[] sigma_mat;
    delete[] epsilon_mat;
    delete[] rc_sqrd_mat;
    delete[] shift_mat;
}


/* ----------------------------------------------------------------------------
   LennardJones destructor, releasing memory of all parameter arrays
------------------------------------------------------------------------------- */

LennardJones::~LennardJones()
{
    free_memory();
}

