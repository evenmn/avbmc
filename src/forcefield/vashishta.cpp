#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

#include "vashishta.h"
#include "../system.h"
#include "../particle.h"


/* ----------------------------------------------------------------------------
   Vashishta constructor, which takes a parameter file 'params'
------------------------------------------------------------------------------- */

Vashishta::Vashishta(System* system_in, const std::string params)
    : ForceField(system_in)
{
    label = "Vashishta";
    paramfile = params;
    read_param_file(params);
    create_label_mapping();
    allocate_memory();
    sort_params();
}


/* ----------------------------------------------------------------------------
   Read parameter file 'params' and store parameters globally. The file takes
   the following form:
   <label1> <label2> <label3> <H> <eta> <Zi> <Zj> <lambda1> <D> <lambda4>
                              <W> <rc> <B> <gamma> <r0> <C> <cos(theta)>
   The file is inspired by LAMMPS, and for an explanation of the 
   symbols and the way it is used, see the LAMMPS Vashishta
   documentation.
------------------------------------------------------------------------------- */

void Vashishta::read_param_file(const std::string params)
{
    nline = 0;
    std::ifstream infile(params);
    if (infile.is_open()) {
        std::string line, label1, label2, label3;
        double H, eta, Zi, Zj, lambda1, D, lambda4, W, rc, B, gamma, r0, C, costheta;
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
            else if (iss >> label1 >> label2 >> label3 >> H >> eta >> Zi >> Zj >> lambda1
                         >> D >> lambda4 >> W >> rc >> B >> gamma >> r0 >> C >> costheta) { 
                label1_vec.push_back(label1);
                label2_vec.push_back(label2);
                label3_vec.push_back(label3);
                H_vec.push_back(H);
                eta_vec.push_back(eta);
                Zi_vec.push_back(Zi);
                Zj_vec.push_back(Zj);
                lambda1_vec.push_back(lambda1);
                D_vec.push_back(D);
                lambda4_vec.push_back(lambda4);
                W_vec.push_back(W);
                rc_vec.push_back(rc);
                B_vec.push_back(B);
                gamma_vec.push_back(gamma);
                r0_vec.push_back(r0);
                C_vec.push_back(C);
                costheta_vec.push_back(costheta);
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
   stored in three-dimensional arrays
------------------------------------------------------------------------------- */

void Vashishta::sort_params()
{
    // link list of chemical symbols to list of type indices
    unsigned int type1, type2, type3, i;
    double lambda1inv, lambda4inv, shift_factor;
    std::vector<int> types1_vec, types2_vec, types3_vec;
    for(std::string label : label1_vec){
        types1_vec.push_back(label2type.at(label));
    }
    for(std::string label : label2_vec){
        types2_vec.push_back(label2type.at(label));
    }
    for(std::string label : label3_vec){
        types3_vec.push_back(label2type.at(label));
    }
    // fill up matrices with parameters
    for (i=0; i < nline; i++) {
        type1 = types1_vec[i];
        type2 = types2_vec[i];
        type3 = types3_vec[i];

        // two-body parameters
        if (type2 == type3) {  
            H_mat[type1][type2] = H_vec[i];
            H_mat[type2][type1] = H_vec[i];
            eta_mat[type1][type2] = eta_vec[i];
            eta_mat[type2][type1] = eta_vec[i];
            Zi_mat[type1][type2] = Zi_vec[i];
            Zi_mat[type2][type1] = Zi_vec[i];
            Zj_mat[type1][type2] = Zj_vec[i];
            Zj_mat[type2][type1] = Zj_vec[i];
            lambda1inv = 1.0 / lambda1_vec[i];
            lambda1inv_mat[type1][type2] = lambda1inv;
            lambda1inv_mat[type2][type1] = lambda1inv;
            D_mat[type1][type2] = D_vec[i];
            D_mat[type2][type1] = D_vec[i];
            lambda4inv = 1.0 / lambda4_vec[i];
            lambda4inv_mat[type1][type2] = lambda4inv;
            lambda4inv_mat[type2][type1] = lambda4inv;
            W_mat[type1][type2] = W_vec[i];
            W_mat[type2][type1] = W_vec[i];
            rc_mat[type1][type2] = rc_vec[i];
            rc_mat[type2][type1] = rc_vec[i];
            gamma_mat[type1][type2] = gamma_vec[i];
            gamma_mat[type2][type1] = gamma_vec[i];
            r0_mat[type1][type2] = r0_vec[i];
            r0_mat[type2][type1] = r0_vec[i];

            // construct shift matrix
            shift_factor = H_vec[i] / std::pow(rc_vec[i], eta_vec[i]);
            shift_factor += Zi_vec[i] * Zj_vec[i] * std::exp(-rc_vec[i] * lambda1inv) / rc_vec[i];
            shift_factor -= D_vec[i] * std::exp(-rc_vec[i] * lambda4inv) / std::pow(rc_vec[i], 4);
            shift_factor -= W_vec[i] / std::pow(rc_vec[i], 6);
            shift_mat[type1][type2] = shift_factor;
            shift_mat[type2][type1] = shift_factor;
        }
        // three-body parameters
        B_mat[type1][type2][type3] = B_vec[i];
        B_mat[type1][type3][type2] = B_vec[i];
        B_mat[type2][type1][type3] = B_vec[i];
        B_mat[type2][type3][type1] = B_vec[i];
        B_mat[type3][type1][type2] = B_vec[i];
        B_mat[type3][type2][type1] = B_vec[i];
        C_mat[type1][type2][type3] = C_vec[i];
        C_mat[type1][type3][type2] = C_vec[i];
        C_mat[type2][type1][type3] = C_vec[i];
        C_mat[type2][type3][type1] = C_vec[i];
        C_mat[type3][type1][type2] = C_vec[i];
        C_mat[type3][type2][type1] = C_vec[i];
        costheta_mat[type1][type2][type3] = costheta_vec[i];
        costheta_mat[type1][type3][type2] = costheta_vec[i];
        costheta_mat[type2][type1][type3] = costheta_vec[i];
        costheta_mat[type2][type3][type1] = costheta_vec[i];
        costheta_mat[type3][type1][type2] = costheta_vec[i];
        costheta_mat[type3][type2][type1] = costheta_vec[i];
    }
}


/* ----------------------------------------------------------------------------
   Compute two-body interaction energy between two particles i and j of types
   'typei' and 'typej', respectively, separated by a distance 'rij'.
------------------------------------------------------------------------------- */

double Vashishta::comp_twobody_par(const int typei, const int typej, const double rij,
                                   std::valarray<double> &force, const bool comp_force)
{
    double rijinv, rijinv2, rijinv4, rijinv6, energy;

    energy = 0.;
    if (rij < rc_mat[typei][typej]) {
        rijinv = 1.0 / rij;
        rijinv2 = rijinv * rijinv;
        rijinv4 = rijinv2 * rijinv2;
        rijinv6 = rijinv4 * rijinv2;
        energy += H_mat[typei][typej] * std::pow(rijinv, eta_mat[typei][typej]);
        energy += Zi_mat[typei][typej] * Zj_mat[typei][typej] * rijinv * std::exp(-rij * lambda1inv_mat[typei][typej]);
        energy -= D_mat[typei][typej] * rijinv4 * std::exp(-rij * lambda4inv_mat[typei][typej]);
        energy -= W_mat[typei][typej] * rijinv6;
        //if (comp_force) {
            // do something
        //}
        energy -= shift_mat[typei][typej];
    }
    return energy;
}


/* ----------------------------------------------------------------------------
   Compute three-body interaction energy between three particles i, j and k of
   types 'typei', 'typej' and 'typek', respectively. 'delij' is the distance
   vector from i to j, while 'delik' is the distance vector from i to k. 'rij'
   is the actual distance between particle i and j.
------------------------------------------------------------------------------- */

double Vashishta::comp_threebody_par(const int typei, const int typej,
        const int typek, const std::valarray<double> delij,
        const std::valarray<double> delik, const double rij,
        std::valarray<double> &force, const bool comp_force)
{
    double rik, expij, expik, costhetaijk, delcos, delcossq, energy;

    energy = 0.;
    rik = std::sqrt(norm(delik));
    if (rij < r0_mat[typei][typej] && rik < r0_mat[typei][typek]) {
        costhetaijk = (delij * delik).sum() / (rij * rik);
        expij = gamma_mat[typei][typej] / (rij - r0_mat[typei][typej]);
        expik = gamma_mat[typei][typek] / (rik - r0_mat[typei][typek]);
        delcos = costheta_mat[typei][typej][typek] - costhetaijk;
        delcossq = delcos * delcos;
        energy += B_mat[typei][typej][typek] * delcossq / 
            (1 + C_mat[typei][typej][typek] * delcossq) * std::exp(expij + expik);
        //if (comp_force) {
        //}
    }
    return (energy);
}


/* ----------------------------------------------------------------------------
   Compute energy contribution of a particle 'i', given a vector containing all
   particles.
------------------------------------------------------------------------------- */

double Vashishta::comp_energy_par(const std::vector<Particle> particles, const int i,
                                  std::valarray<double> &force, const bool comp_force)
{
    // declare variables
    double rij, energy;
    int typei, typej, typek, npar, j, k;
    std::valarray<double> delij, delik;
    typei = particles[i].type; 
    npar = particles.size();

    energy = 0.;
    for(j=0; j<i; j++){
        // two-body
        typej = particles[j].type;
        delij = particles[j].r - particles[i].r;
        rij = std::sqrt(norm(delij));
        energy += comp_twobody_par(typei, typej, rij, force, comp_force);

        for(k=0; k<npar; k++){
            if (k==i || k == j) continue;

            // three-body
            typek = particles[k].type;
            delik = particles[k].r - particles[i].r;
            energy += comp_threebody_par(typei, typej, typek, delij, delik, rij, force, comp_force);
        }
    }

    for(j=i+1; j<npar; j++){
        // two-body
        typej = particles[j].type;
        delij = particles[j].r - particles[i].r;
        rij = std::sqrt(norm(delij));
        energy += comp_twobody_par(typei, typej, rij, force, comp_force);

        for(k=0; k<npar; k++){
            if (k==i || k == j) continue;

            // three-body
            typek = particles[k].type;
            delik = particles[k].r - particles[i].r;
            energy += comp_threebody_par(typei, typej, typek, delij, delik, rij, force, comp_force);
        }
    }
    return energy;
}


/* ----------------------------------------------------------------------------
   Compute total energy of a system consisting of a set of particles stored
   in the vector 'particles'. This is just used to verify that the energy 
   difference is computed correctly.
------------------------------------------------------------------------------- */
/*
double Vashishta::comp_energy_tot(const std::vector<Particle> particles)
{
    // declare variables
    double rij, energy;
    int typei, typej, typek, npar, i, j, k;
    std::valarray<double> delij, delik;
    npar = particles.size();

    energy = 0.;
    for (i=0; i<npar; i++) {
        typei = particles[i].type;
        for (j=0; j<i; j++) {
            // two-body
            typej = particles[j].type;
            delij = particles[j].r - particles[i].r;
            rij = std::sqrt(norm(delij));
            energy += comp_twobody_par(typei, typej, rij, force, comp_force);

            for (k=0; k<npar; k++) {
                if (k==i || k==j) continue;

                // three-body
                typek = particles[k].type;
                delik = particles[k].r - particles[i].r;
                energy += comp_threebody_par(typei, typej, typek, delij, delik, rij, force, comp_force);
            }
        }

        for (j=i+1; j<npar; j++) {
            for (k=0; k<npar; k++) {
                if (k==i || k==j) continue;

                // three-body
                typek = particles[k].type;
                delik = particles[k].r - particles[i].r;
                energy += comp_threebody_par(typei, typej, typek, delij, delik, rij, force, comp_force);
            }
        }
    }
    return energy;
}
*/

/* ----------------------------------------------------------------------------
   Allocate memory
------------------------------------------------------------------------------- */

void Vashishta::allocate_memory()
{
    H_mat = new double*[ntype];
    eta_mat = new double*[ntype];
    Zi_mat = new double*[ntype];
    Zj_mat = new double*[ntype];
    lambda1inv_mat = new double*[ntype];
    D_mat = new double*[ntype];
    lambda4inv_mat = new double*[ntype];
    W_mat = new double*[ntype];
    rc_mat = new double*[ntype];
    gamma_mat = new double*[ntype];
    r0_mat = new double*[ntype];
    shift_mat = new double*[ntype];
    B_mat = new double**[ntype];
    C_mat = new double**[ntype];
    costheta_mat = new double**[ntype];
    for (unsigned int i=0; i<ntype; i++) {
        H_mat[i] = new double[ntype];
        eta_mat[i] = new double[ntype];
        Zi_mat[i] = new double[ntype];
        Zj_mat[i] = new double[ntype];
        lambda1inv_mat[i] = new double[ntype];
        D_mat[i] = new double[ntype];
        lambda4inv_mat[i] = new double[ntype];
        W_mat[i] = new double[ntype];
        rc_mat[i] = new double[ntype];
        gamma_mat[i] = new double[ntype];
        r0_mat[i] = new double[ntype];
        shift_mat[i] = new double[ntype];
        B_mat[i] = new double*[ntype];
        C_mat[i] = new double*[ntype];
        costheta_mat[i] = new double*[ntype];
        for (unsigned int j=0; j<ntype; j++) {
            B_mat[i][j] = new double[ntype];
            C_mat[i][j] = new double[ntype];
            costheta_mat[i][j] = new double[ntype];
        }
    }
}


/* ----------------------------------------------------------------------------
   Free memory
------------------------------------------------------------------------------- */

void Vashishta::free_memory()
{
    for (unsigned int i = 0; i < ntype; i++) {
        for (unsigned int j = 0; j < ntype; j++) {
            delete[] B_mat[i][j];
            delete[] C_mat[i][j];
            delete[] costheta_mat[i][j];
        }
        delete[] H_mat[i];
        delete[] eta_mat[i];
        delete[] Zi_mat[i];
        delete[] Zj_mat[i];
        delete[] lambda1inv_mat[i];
        delete[] D_mat[i];
        delete[] lambda4inv_mat[i];
        delete[] W_mat[i];
        delete[] rc_mat[i];
        delete[] B_mat[i];
        delete[] gamma_mat[i];
        delete[] r0_mat[i];
        delete[] shift_mat[i];
        delete[] C_mat[i];
        delete[] costheta_mat[i];
    }
    delete[] H_mat;
    delete[] eta_mat;
    delete[] Zi_mat;
    delete[] Zj_mat;
    delete[] lambda1inv_mat;
    delete[] D_mat;
    delete[] lambda4inv_mat;
    delete[] W_mat;
    delete[] rc_mat;
    delete[] B_mat;
    delete[] gamma_mat;
    delete[] r0_mat;
    delete[] shift_mat;
    delete[] C_mat;
    delete[] costheta_mat;
}


/* ----------------------------------------------------------------------------
   Vashishta destructor, memory of all parameter arrays
   is released.
------------------------------------------------------------------------------- */

Vashishta::~Vashishta()
{
    free_memory();
}

