#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cassert>

#include "lennardjones.h"
#include "../system.h"
#include "../particle.h"


/* -----------------------------------------------------
   This is the default constructor when parameter
   file is not given. It only works for pure Argon
-------------------------------------------------------- */
/*
LennardJones::LennardJones(System* system_in)
    : ForceField(system_in)
{
    label1_vec = {"Ar"};
    label2_vec = {"Ar"};
    sigma_vec = {1.};
    epsilon_vec = {1.};
    rc_vec = {5.};
    nline = 1;
    label = "Lennard-Jones";
    paramfile = "";
}
*/

/* ------------------------------------------------------
   This is the default constructor when a parameter
   file 'params' is given.
--------------------------------------------------------- */

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


/* ------------------------------------------------------
   Read parameter file 'params' and store parameters 
   globally. It takes the following form:
       <label1> <label2> <sigma> <epsilon> <rc>
--------------------------------------------------------- */

void LennardJones::read_param_file(const std::string params)
{
    nline = 0;
    std::ifstream infile(params);
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


/* ---------
   Allocate memory for matrices
---------------------------------- */

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


/* --------------
   Free memory for matrices
------------------------------- */

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

/* ------------------------------------------------------
   Parameters have to be sorted with respect to
   the particle types, and are stored in matrices.

   TODO: Since the possible particle labels (chemical
   elements) are defined by the parameter file, this 
   function can be called by the constructor.
--------------------------------------------------------- */

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


/* ------------------------------------------------------
   Compute interaction energy between two particles of
   types 'typei' and 'typej', respectively, separated
   by a distance vector 'delij'. Updates a force array
   'force' if 'comp_force' is true.
--------------------------------------------------------- */

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


/* -------------------------------------------------------
   Compute energy contribution from a particle 'i'
---------------------------------------------------------- */

double LennardJones::comp_energy_par(const std::vector<Particle> particles, const int i,
                                     std::valarray<double> &force, const bool comp_force)
{
    // declare variables
    int npar, typei, typej, j;
    npar = particles.size();
    typei = particles[i].type;
    std::valarray<double> delij;
    force.resize(system->ndim, 0.);

    double energy = 0.;
    for (j=0; j<i; j++) {
        typej = particles[j].type;
        delij = particles[j].r - particles[i].r;
        energy += comp_twobody_par(typei, typej, delij, force, comp_force);
    }
    for (j=i+1; j<npar; j++) {
        typej = particles[j].type;
        delij = particles[j].r - particles[i].r;
        energy += comp_twobody_par(typei, typej, delij, force, comp_force);
    }
    force *= 24;
    return (4 * energy);
}


/* --------------------------------------------------
   LennardJones destructor, releasing memory of all
   parameter arrays
----------------------------------------------------- */

LennardJones::~LennardJones()
{
    free_memory();
}

