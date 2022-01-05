#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <cassert>
#include <memory>

#include "vashishta.h"
#include "../system.h"
#include "../particle.h"
#include "../molecule.h"


/* -----------------------------------------------------------------
   Vashishta constructor, which takes a parameter file 'params'
-------------------------------------------------------------------- */

Vashishta::Vashishta(System* system_in, const std::string params)
    : ForceField(system_in)
{
    read_param_file(params);
    label = "Vashishta";
    paramfile = params;
}


/* -----------------------------------------------------------------
   Read parameter file 'params' and store parameters
   globally. The file takes the following form:
   <label1> <label2> <label3> <H> <eta> <Zi> <Zj> <lambda1> <D> <lambda4>
                              <W> <rc> <B> <gamma> <r0> <C> <cos(theta)>
   The file is inspired by LAMMPS, and for an explanation of the 
   symbols and the way it is used, see the LAMMPS Vashishta
   documentation.
-------------------------------------------------------------------- */

void Vashishta::read_param_file(const std::string params)
{
    nline = 0;
    std::ifstream infile(params);
    std::string line, label1, label2, label3;
    double H, eta, Zi, Zj, lambda1, D, lambda4, W, rc, B, gamma, r0, C, costheta;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        if (line.rfind("#", 0) == 0){
            // comments are allowed in parameter file
            continue;
        }
        else if (line.empty()){
            // empty lines are allowed in parameter file
            continue;
        }
        else if (iss >> label1 >> label2 >> label3 >> H >> eta >> Zi >> Zj >> lambda1
                     >> D >> lambda4 >> W >> rc >> B >> gamma >> r0 >> C >> costheta){ 
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
        else{
            std::cout << "Warning: Corrupt line in parameter file!" << std::endl;
            std::cout << "Ignoring line: '" + line + "'" << std::endl;
        }
    }
}


/* -----------------------------------------------------------------
   Parameters have to be sorted with respect to the particle
   types, and are stored in three-dimensional arrays
-------------------------------------------------------------------- */

void Vashishta::sort_params()
{
    // link list of chemical symbols to list of type indices
    std::vector<int> types1_vec;
    std::vector<int> types2_vec;
    std::vector<int> types3_vec;
    for(std::string label : label1_vec){
        types1_vec.push_back(system->label2type.at(label));
    }
    for(std::string label : label2_vec){
        types2_vec.push_back(system->label2type.at(label));
    }
    for(std::string label : label3_vec){
        types3_vec.push_back(system->label2type.at(label));
    }
    /*
    for(std::string label : label1_vec){
        bool assigned = false;
        for(int j=0; j<system->ntype; j++){
            if(label == system->unique_labels[j]){
                types1_vec.push_back(j);
                assigned = true;
            }
        }
        assert(assigned);
    }
    for(std::string label : label1_vec){
        types2_vec.push_back(system->label2type.at(label));
    }
    for(std::string label : label2_vec){
        bool assigned = false;
        for(int j=0; j<system->ntype; j++){
            if(label == system->unique_labels[j]){
                types2_vec.push_back(j);
                assigned = true;
            }
        }
        assert(assigned);
    }
    for(std::string label : label1_vec){
        types1_vec.push_back(system->label2type.at(label));
    }

    for(std::string label : label3_vec){
        bool assigned = false;
        for(int j=0; j<system->ntype; j++){
            if(label == system->unique_labels[j]){
                types3_vec.push_back(j);
                assigned = true;
            }
        }
        assert(assigned);
    }
    */
    // allocate memory for matrices
    H_mat = new double*[system->ntype];
    eta_mat = new double*[system->ntype];
    Zi_mat = new double*[system->ntype];
    Zj_mat = new double*[system->ntype];
    lambda1inv_mat = new double*[system->ntype];
    D_mat = new double*[system->ntype];
    lambda4inv_mat = new double*[system->ntype];
    W_mat = new double*[system->ntype];
    rc_mat = new double*[system->ntype];
    gamma_mat = new double*[system->ntype];
    r0_mat = new double*[system->ntype];
    B_mat = new double**[system->ntype];
    C_mat = new double**[system->ntype];
    costheta_mat = new double**[system->ntype];
    for(int i=0; i<system->ntype; i++){
        H_mat[i] = new double[system->ntype];
        eta_mat[i] = new double[system->ntype];
        Zi_mat[i] = new double[system->ntype];
        Zj_mat[i] = new double[system->ntype];
        lambda1inv_mat[i] = new double[system->ntype];
        D_mat[i] = new double[system->ntype];
        lambda4inv_mat[i] = new double[system->ntype];
        W_mat[i] = new double[system->ntype];
        rc_mat[i] = new double[system->ntype];
        gamma_mat[i] = new double[system->ntype];
        r0_mat[i] = new double[system->ntype];
        B_mat[i] = new double*[system->ntype];
        C_mat[i] = new double*[system->ntype];
        costheta_mat[i] = new double*[system->ntype];
        for(int j=0; j<system->ntype; j++){
            B_mat[i][j] = new double[system->ntype];
            C_mat[i][j] = new double[system->ntype];
            costheta_mat[i][j] = new double[system->ntype];
        }
    }
    // fill up matrices with parameters
    int type1, type2, type3;
    for(int i=0; i<nline; i++){
        type1 = types1_vec[i];
        type2 = types2_vec[i];
        type3 = types3_vec[i];

        // two-body parameters
        if (type2 == type3){  
            H_mat[type1][type2] = H_vec[i];
            H_mat[type2][type1] = H_vec[i];
            eta_mat[type1][type2] = eta_vec[i];
            eta_mat[type2][type1] = eta_vec[i];
            Zi_mat[type1][type2] = Zi_vec[i];
            Zi_mat[type2][type1] = Zi_vec[i];
            Zj_mat[type1][type2] = Zj_vec[i];
            Zj_mat[type2][type1] = Zj_vec[i];
            lambda1inv_mat[type1][type2] = 1.0 / lambda1_vec[i];
            lambda1inv_mat[type2][type1] = 1.0 / lambda1_vec[i];
            D_mat[type1][type2] = D_vec[i];
            D_mat[type2][type1] = D_vec[i];
            lambda4inv_mat[type1][type2] = 1.0 / lambda4_vec[i];
            lambda4inv_mat[type2][type1] = 1.0 / lambda4_vec[i];
            W_mat[type1][type2] = W_vec[i];
            W_mat[type2][type1] = W_vec[i];
            rc_mat[type1][type2] = rc_vec[i];
            rc_mat[type2][type1] = rc_vec[i];
            gamma_mat[type1][type2] = gamma_vec[i];
            gamma_mat[type2][type1] = gamma_vec[i];
            r0_mat[type1][type2] = r0_vec[i];
            r0_mat[type2][type1] = r0_vec[i];
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


/* ----------------------------------------------------------------
   Compute two-body interaction energy between two particles
   i and j of types 'typei' and 'typej', respectively, 
   separated by a distance 'rij'.
------------------------------------------------------------------- */

double Vashishta::comp_twobody_par(const int typei, const int typej, const double rij)
{
    double rijinv, rijinv2, rijinv4, rijinv6, energy;

    energy = 0.0;
    if(rij < rc_mat[typei][typej]){
        rijinv = 1.0 / rij;
        rijinv2 = rijinv * rijinv;
        rijinv4 = rijinv2 * rijinv2;
        rijinv6 = rijinv4 * rijinv2;
        energy += H_mat[typei][typej] * std::pow(rijinv, eta_mat[typei][typej]);
        energy += Zi_mat[typei][typej] * Zj_mat[typei][typej] * rijinv * std::exp(-rij * lambda1inv_mat[typei][typej]);
        energy -= D_mat[typei][typej] * rijinv4 * std::exp(-rij * lambda4inv_mat[typei][typej]);
        energy -= W_mat[typei][typej] * rijinv6;
    }
    return energy;
}


/* ---------------------------------------------------------------
   Compute three-body interaction energy between three 
   particles i, j and k of types 'typei', 'typej' and 'typek',
   respectively. 'delij' is the distance vector from i to j,
   while 'delik' is the distance vector from i to k. 'rij'
   is the actual distance between particle i and j.
------------------------------------------------------------------ */

double Vashishta::comp_threebody_par(const int typei, const int typej, const int typek,
        const std::valarray<double> delij, const std::valarray<double> delik, const double rij)
{
    double rik, expij, expik, costhetaijk, delcos, delcossq, energy;

    rik = std::sqrt(norm(delik));

    energy = 0.0;
    if(rij < r0_mat[typei][typej] && rik < r0_mat[typei][typek]){
        costhetaijk = (delij * delik).sum() / (rij * rik);
        expij = std::exp(gamma_mat[typei][typej] / (rij - r0_mat[typei][typej]));
        expik = std::exp(gamma_mat[typei][typek] / (rik - r0_mat[typei][typek]));
        delcos = costheta_mat[typei][typej][typek] - costhetaijk;
        delcossq = delcos * delcos;
        energy += B_mat[typei][typej][typek] * delcossq / (1 + C_mat[typei][typej][typek] * delcossq) * expij * expik;
    }
    return energy;
}


/* ---------------------------------------------------------------
   Compute energy contribution from molecule consisting of
   one or several atoms. This function might be redundant.
------------------------------------------------------------------ */

/*
double Vashishta::comp_energy_mol(const std::vector<Particle *> particles, Molecule* molecule)
{
    double energy = 0.0;
    for (int i : molecule->atoms_idx){
        energy += comp_energy_par(particles, i);
    }
    return energy;
}
*/
double Vashishta::comp_energy_mol(const std::vector<std::shared_ptr<Particle> > particles, Molecule* molecule)
{
    double energy = 0.0;
    for (int i : molecule->atoms_idx){
        energy += comp_energy_par(particles, i);
    }
    return energy;
}


/* ---------------------------------------------------------------
   Compute energy contribution of a particle 'i', given a 
   vector containing all particles.
------------------------------------------------------------------ */
/*
double Vashishta::comp_energy_par(const std::vector<Particle *> particles, const int i)
{
    // declare variables
    int typei = particles[i]->type; 
    int typej, typek;
    int npar = particles.size();
    double rij, energy;
    std::valarray<double> delij, delik;

    energy = 0.0;
    for(int j=0; j<i; j++){
        // two-body
        typej = particles[j]->type;
        delij = particles[j]->r - particles[i]->r;
        rij = std::sqrt(norm(delij));
        energy += comp_twobody_par(typei, typej, rij);

        for(int k=0; k<npar; k++){
            if (k==i || k == j) continue;

            // three-body
            typek = particles[k]->type;
            delik = particles[k]->r - particles[i]->r;
            energy += comp_threebody_par(typei, typej, typek, delij, delik, rij);
        }
    }

    for(int j=i+1; j<npar; j++){
        // two-body
        typej = particles[j]->type;
        delij = particles[j]->r - particles[i]->r;
        rij = std::sqrt(norm(delij));  // might be a good idea to write this out instead
        energy += comp_twobody_par(typei, typej, rij);

        for(int k=0; k<npar; k++){
            if (k==i || k == j) continue;

            // three-body
            typek = particles[k]->type;
            delik = particles[k]->r - particles[i]->r;
            energy += comp_threebody_par(typei, typej, typek, delij, delik, rij);
        }
    }
    return (energy);
}
*/
double Vashishta::comp_energy_par(const std::vector<std::shared_ptr<Particle> > particles, const int i)
{
    // declare variables
    int typei = particles[i]->type; 
    int typej, typek;
    int npar = particles.size();
    double rij, energy;
    std::valarray<double> delij, delik;

    energy = 0.0;
    for(int j=0; j<i; j++){
        // two-body
        typej = particles[j]->type;
        delij = particles[j]->r - particles[i]->r;
        rij = std::sqrt(norm(delij));
        energy += comp_twobody_par(typei, typej, rij);

        for(int k=0; k<npar; k++){
            if (k==i || k == j) continue;

            // three-body
            typek = particles[k]->type;
            delik = particles[k]->r - particles[i]->r;
            energy += comp_threebody_par(typei, typej, typek, delij, delik, rij);
        }
    }

    for(int j=i+1; j<npar; j++){
        // two-body
        typej = particles[j]->type;
        delij = particles[j]->r - particles[i]->r;
        rij = std::sqrt(norm(delij));  // might be a good idea to write this out instead
        energy += comp_twobody_par(typei, typej, rij);

        for(int k=0; k<npar; k++){
            if (k==i || k == j) continue;

            // three-body
            typek = particles[k]->type;
            delik = particles[k]->r - particles[i]->r;
            energy += comp_threebody_par(typei, typej, typek, delij, delik, rij);
        }
    }
    return (energy);
}

/* ---------------------------------------------------------
   Vashishta destructor, memory of all parameter arrays
   is released. Do this in a more elegant way in the 
   future
------------------------------------------------------------ */

Vashishta::~Vashishta()
{
    for(int i = 0; i < system->ntype; i++){
        for(int j = 0; j < system->ntype; j++){
            delete B_mat[i][j];
            delete C_mat[i][j];
            delete costheta_mat[i][j];
        }
        delete H_mat[i];
        delete eta_mat[i];
        delete Zi_mat[i];
        delete Zj_mat[i];
        delete lambda1inv_mat[i];
        delete D_mat[i];
        delete lambda4inv_mat[i];
        delete W_mat[i];
        delete rc_mat[i];
        delete B_mat[i];
        delete gamma_mat[i];
        delete r0_mat[i];
        delete C_mat[i];
        delete costheta_mat[i];
    }
    delete H_mat;
    delete eta_mat;
    delete Zi_mat;
    delete Zj_mat;
    delete lambda1inv_mat;
    delete D_mat;
    delete lambda4inv_mat;
    delete W_mat;
    delete rc_mat;
    delete B_mat;
    delete gamma_mat;
    delete r0_mat;
    delete C_mat;
    delete costheta_mat;
}

