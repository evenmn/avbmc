#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "vashishta.h"
#include "../box.h"
#include "../particle.h"


/* -----------------------------------------------------------------
   Vashishta constructor, which takes a parameter file 'params'
-------------------------------------------------------------------- */

Vashishta::Vashishta(Box* box_in, const std::string params)
    : ForceField(box_in)
{
    read_param_file(params);
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
            std::cout << "Warning: Corrupt line in parameter file" << std::endl;
            std::cout << "Line: " + line << std::endl;
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
        bool assigned = false;
        for(int j=0; j<box->ntype; j++){
            if(label == box->unique_labels[j]){
                types1_vec.push_back(j);
                assigned = true;
            }
        }
        assert(assigned);
    }
    for(std::string label : label2_vec){
        bool assigned = false;
        for(int j=0; j<box->ntype; j++){
            if(label == box->unique_labels[j]){
                types2_vec.push_back(j);
                assigned = true;
            }
        }
        assert(assigned);
    }

    for(std::string label : label3_vec){
        bool assigned = false;
        for(int j=0; j<box->ntype; j++){
            if(label == box->unique_labels[j]){
                types3_vec.push_back(j);
                assigned = true;
            }
        }
        assert(assigned);
    }

    // allocate memory for matrices
    H_mat = new double**[box->ntype];
    eta_mat = new double**[box->ntype];
    Zi_mat = new double**[box->ntype];
    Zj_mat = new double**[box->ntype];
    lambda1_mat = new double**[box->ntype];
    D_mat = new double**[box->ntype];
    lambda4_mat = new double**[box->ntype];
    W_mat = new double**[box->ntype];
    rc_sqrd_mat = new double**[box->ntype];
    B_mat = new double**[box->ntype];
    gamma_mat = new double**[box->ntype];
    r0_mat = new double**[box->ntype];
    C_mat = new double**[box->ntype];
    costheta_mat = new double**[box->ntype];
    for(int i=0; i<box->ntype; i++){
        H_mat[i] = new double*[box->ntype];
        eta_mat[i] = new double*[box->ntype];
        Zi_mat[i] = new double*[box->ntype];
        Zj_mat[i] = new double*[box->ntype];
        lambda1_mat[i] = new double*[box->ntype];
        D_mat[i] = new double*[box->ntype];
        lambda4_mat[i] = new double*[box->ntype];
        W_mat[i] = new double*[box->ntype];
        rc_sqrd_mat[i] = new double*[box->ntype];
        B_mat[i] = new double*[box->ntype];
        gamma_mat[i] = new double*[box->ntype];
        r0_mat[i] = new double*[box->ntype];
        C_mat[i] = new double*[box->ntype];
        costheta_mat[i] = new double*[box->ntype];
        for(int j=0; j<box->ntype; j++){
            H_mat[i][j] = new double[box->ntype];
            eta_mat[i][j] = new double[box->ntype];
            Zi_mat[i][j] = new double[box->ntype];
            Zj_mat[i][j] = new double[box->ntype];
            lambda1_mat[i][j] = new double[box->ntype];
            D_mat[i][j] = new double[box->ntype];
            lambda4_mat[i][j] = new double[box->ntype];
            W_mat[i][j] = new double[box->ntype];
            rc_sqrd_mat[i][j] = new double[box->ntype];
            B_mat[i][j] = new double[box->ntype];
            gamma_mat[i][j] = new double[box->ntype];
            r0_mat[i][j] = new double[box->ntype];
            C_mat[i][j] = new double[box->ntype];
            costheta_mat[i][j] = new double[box->ntype];
        }
    }
    // fill up matrices with parameters
    int type1, type2, type3;
    for(int i=0; i<nline; i++){
        type1 = types1_vec[i];
        type2 = types2_vec[i];
        type3 = types3_vec[i];
        H_mat[type1][type2][type3] = H_vec[i];
        H_mat[type1][type3][type2] = H_vec[i];
        H_mat[type2][type1][type3] = H_vec[i];
        H_mat[type2][type3][type1] = H_vec[i];
        H_mat[type3][type1][type2] = H_vec[i];
        H_mat[type3][type2][type1] = H_vec[i];
        eta_mat[type1][type2][type3] = eta_vec[i];
        eta_mat[type1][type3][type2] = eta_vec[i];
        eta_mat[type2][type1][type3] = eta_vec[i];
        eta_mat[type2][type3][type1] = eta_vec[i];
        eta_mat[type3][type1][type2] = eta_vec[i];
        eta_mat[type3][type2][type1] = eta_vec[i];
        Zi_mat[type1][type2][type3] = Zi_vec[i];
        Zi_mat[type1][type3][type2] = Zi_vec[i];
        Zi_mat[type2][type1][type3] = Zi_vec[i];
        Zi_mat[type2][type3][type1] = Zi_vec[i];
        Zi_mat[type3][type1][type2] = Zi_vec[i];
        Zi_mat[type3][type2][type1] = Zi_vec[i];
        Zj_mat[type1][type2][type3] = Zj_vec[i];
        Zj_mat[type1][type3][type2] = Zj_vec[i];
        Zj_mat[type2][type1][type3] = Zj_vec[i];
        Zj_mat[type2][type3][type1] = Zj_vec[i];
        Zj_mat[type3][type1][type2] = Zj_vec[i];
        Zj_mat[type3][type2][type1] = Zj_vec[i];
        lambda1_mat[type1][type2][type3] = lambda1_vec[i];
        lambda1_mat[type1][type3][type2] = lambda1_vec[i];
        lambda1_mat[type2][type1][type3] = lambda1_vec[i];
        lambda1_mat[type2][type3][type1] = lambda1_vec[i];
        lambda1_mat[type3][type1][type2] = lambda1_vec[i];
        lambda1_mat[type3][type2][type1] = lambda1_vec[i];
        D_mat[type1][type2][type3] = D_vec[i];
        D_mat[type1][type3][type2] = D_vec[i];
        D_mat[type2][type1][type3] = D_vec[i];
        D_mat[type2][type3][type1] = D_vec[i];
        D_mat[type3][type1][type2] = D_vec[i];
        D_mat[type3][type2][type1] = D_vec[i];
        lambda4_mat[type1][type2][type3] = lambda4_vec[i];
        lambda4_mat[type1][type3][type2] = lambda4_vec[i];
        lambda4_mat[type2][type1][type3] = lambda4_vec[i];
        lambda4_mat[type2][type3][type1] = lambda4_vec[i];
        lambda4_mat[type3][type1][type2] = lambda4_vec[i];
        lambda4_mat[type3][type2][type1] = lambda4_vec[i];
        W_mat[type1][type2][type3] = W_vec[i];
        W_mat[type1][type3][type2] = W_vec[i];
        W_mat[type2][type1][type3] = W_vec[i];
        W_mat[type2][type3][type1] = W_vec[i];
        W_mat[type3][type1][type2] = W_vec[i];
        W_mat[type3][type2][type1] = W_vec[i];
        double rc_sqrd = rc_vec[i] * rc_vec[i];
        rc_sqrd_mat[type1][type2][type3] = rc_sqrd;
        rc_sqrd_mat[type1][type3][type2] = rc_sqrd;
        rc_sqrd_mat[type2][type1][type3] = rc_sqrd;
        rc_sqrd_mat[type2][type3][type1] = rc_sqrd;
        rc_sqrd_mat[type3][type1][type2] = rc_sqrd;
        rc_sqrd_mat[type3][type2][type1] = rc_sqrd;
        B_mat[type1][type2][type3] = B_vec[i];
        B_mat[type1][type3][type2] = B_vec[i];
        B_mat[type2][type1][type3] = B_vec[i];
        B_mat[type2][type3][type1] = B_vec[i];
        B_mat[type3][type1][type2] = B_vec[i];
        B_mat[type3][type2][type1] = B_vec[i];
        gamma_mat[type1][type2][type3] = gamma_vec[i];
        gamma_mat[type1][type3][type2] = gamma_vec[i];
        gamma_mat[type2][type1][type3] = gamma_vec[i];
        gamma_mat[type2][type3][type1] = gamma_vec[i];
        gamma_mat[type3][type1][type2] = gamma_vec[i];
        gamma_mat[type3][type2][type1] = gamma_vec[i];
        r0_mat[type1][type2][type3] = r0_vec[i];
        r0_mat[type1][type3][type2] = r0_vec[i];
        r0_mat[type2][type1][type3] = r0_vec[i];
        r0_mat[type2][type3][type1] = r0_vec[i];
        r0_mat[type3][type1][type2] = r0_vec[i];
        r0_mat[type3][type2][type1] = r0_vec[i];
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
   Compute two-body contribution to total energy for a particle 'i'
------------------------------------------------------------------- */

double Vashishta::comp_twobody_par(const std::vector<Particle *> particles, const int i)
{
    int typei = particles[i]->type;
    int typej;
    double sdist;

    double energy = 0;
    for(int j=0; j<box->npar; j++){
        sdist = std::pow(particles[j]->r - particles[i]->r, 2).sum();
        typej = particles[j]->type;
        energy += 1;
    }
}


/* ---------------------------------------------------------------
   Compute energy contribution from molecule consisting of
   one or several atoms
------------------------------------------------------------------ */
double Vashishta::comp_energy_mol(std::vector<Particle *> particles, Molecule* molecule)
{
    return 0.;
}


/* ---------------------------------------------------------------
   Compute energy contribution of a particle i
------------------------------------------------------------------ */

double Vashishta::comp_energy_par(const std::vector<Particle *> particles, const int i)
{
    // declare variables
    int typei = particles[i]->type; 
    int typej, typek;
    double rij, rik, rijsq, riksq, expij, expik, rijinv, rijinv4, rijinv6;
    double costhetaijk, delcos, delcossq;
    std::valarray<double> delij, delik;

    double energy = 0;
    for(int j=0; j<i; j++){

        // temporary ij quantities
        typej = particles[j]->type;
        delij = particles[j]->r - particles[i]->r;
        rijsq = std::pow(delij, 2).sum();
        rij = std::sqrt(rijsq);
        rijinv = 1.0 / rij;
        rijinv4 = 1.0 / (rijsq * rijsq);
        rijinv6 = rijinv4 / rijsq;
        expij = std::exp(gamma_mat[typei][typej][typej] / (rij - r0_mat[typei][typej][typej]));

        if(rijsq < rc_sqrd_mat[typei][typej][typej]){
            energy += H_mat[typei][typej][typej] * std::pow(rijinv, eta_mat[typei][typej][typej]);
            energy += Zi_mat[typei][typej][typej] * Zj_mat[typei][typej][typej] * rijinv * std::exp(-rij/lambda1_mat[typei][typej][typej]);
            energy -= D_mat[typei][typej][typej] * rijinv4 * std::exp(-rij/lambda4_mat[typei][typej][typek]);
            energy -= W_mat[typei][typej][typej] * rijinv6;
        }
        for(int k=0; k<box->npar; k++){

            // temporary ik quantities
            typek = particles[k]->type;
            delik = particles[k]->r - particles[i]->r;
            riksq = std::pow(delik, 2).sum();
            rik = std::sqrt(riksq);
            expik = std::exp(gamma_mat[typei][typek][typek] / (rik - r0_mat[typei][typek][typek]));

            if(rij < r0_mat[typei][typej][typej] && rik < r0_mat[typei][typek][typek]){
                costhetaijk = (delij * delik).sum() / (rij * rik);
                delcos = costheta_mat[typei][typej][typek] - costhetaijk;
                delcossq = delcos * delcos;
                energy += B_mat[typei][typej][typek] * delcossq / (1 + C_mat[typei][typej][typek] * delcossq) * expij * expik;
            }
        }
    }
    for(int j=i+1; j<box->npar; j++){

        // temporary ij quantities
        typej = particles[j]->type;
        delij = particles[j]->r - particles[i]->r;
        rijsq = std::pow(delij, 2).sum();
        rij = std::sqrt(rijsq);
        rijinv = 1.0 / rij;
        rijinv4 = 1.0 / (rijsq * rijsq);
        rijinv6 = rijinv4 / rijsq;
        expij = std::exp(gamma_mat[typei][typej][typej] / (rij - r0_mat[typei][typej][typej]));

        if(rijsq < rc_sqrd_mat[typei][typej][typej]){
            energy += H_mat[typei][typej][typej] * std::pow(rijinv, eta_mat[typei][typej][typej]);
            energy += Zi_mat[typei][typej][typej] * Zj_mat[typei][typej][typej] * rijinv * std::exp(-rij/lambda1_mat[typei][typej][typej]);
            energy -= D_mat[typei][typej][typej] * rijinv4 * std::exp(-rij/lambda4_mat[typei][typej][typek]);
            energy -= W_mat[typei][typej][typej] * rijinv6;
        }
        for(int k=0; k<box->npar; k++){

            // temporary ik quantities
            typek = particles[k]->type;
            delik = particles[k]->r - particles[i]->r;
            riksq = std::pow(delik, 2).sum();
            rik = std::sqrt(riksq);
            expik = std::exp(gamma_mat[typei][typek][typek] / (rik - r0_mat[typei][typek][typek]));

            if(rij < r0_mat[typei][typej][typej] && rik < r0_mat[typei][typek][typek]){
                costhetaijk = (delij * delik).sum() / (rij * rik);
                delcos = costheta_mat[typei][typej][typek] - costhetaijk;
                delcossq = delcos * delcos;
                energy += B_mat[typei][typej][typek] * delcossq / (1 + C_mat[typei][typej][typek] * delcossq) * expij * expik;
            }
        }
    }
    return (energy);
}


/* ---------------------------------------------------------
   Vashishta destructor, memory of all parameter arrays
   are released. Do this in a more elegant way in the 
   future
------------------------------------------------------------ */

Vashishta::~Vashishta()
{
    for(int i = 0; i < box->ntype; i++){
        for(int j = 0; j < box->ntype; j++){
            delete H_mat[i][j];
            delete eta_mat[i][j];
            delete Zi_mat[i][j];
            delete Zj_mat[i][j];
            delete lambda1_mat[i][j];
            delete D_mat[i][j];
            delete lambda4_mat[i][j];
            delete W_mat[i][j];
            delete rc_sqrd_mat[i][j];
            delete B_mat[i][j];
            delete gamma_mat[i][j];
            delete r0_mat[i][j];
            delete C_mat[i][j];
            delete costheta_mat[i][j];
        }
        delete H_mat[i];
        delete eta_mat[i];
        delete Zi_mat[i];
        delete Zj_mat[i];
        delete lambda1_mat[i];
        delete D_mat[i];
        delete lambda4_mat[i];
        delete W_mat[i];
        delete rc_sqrd_mat[i];
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
    delete lambda1_mat;
    delete D_mat;
    delete lambda4_mat;
    delete W_mat;
    delete rc_sqrd_mat;
    delete B_mat;
    delete gamma_mat;
    delete r0_mat;
    delete C_mat;
    delete costheta_mat;
}

