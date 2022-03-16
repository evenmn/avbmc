#include <iostream>
#include <string>
#include <vector>
#include <valarray>

#include "forcefield.h"
#include "../box.h"
#include "../system.h"
#include "../particle.h"

ForceField::ForceField(Box* box_in)
{
    box = box_in;
    ntype = nline = 0;

    if (!box->store_distance) {
        comp_energy_par = &ForceField::comp_energy_par_neigh0_eng0;
    }
    else if (!box->store_energy) {
        comp_energy_par = &ForceField::comp_energy_par_neigh1_eng0;
    }
    else {
        comp_energy_par = &ForceField::comp_energy_par_neigh1_eng1;
    }
}


/* -------------------------------------------------------
   Compute energy contribution from a particle 'i'
---------------------------------------------------------- */

double ForceField::comp_energy_par_force0(const int i)
{
    std::valarray<double> force;
    return (this->*(this->comp_energy_par))(i, force, false);
}


/* ---------------------
   Initialize energy and force matrices
---- */

void ForceField::initialize()
{
    unsigned int npar, i, j;

    npar = box->npar;
    poteng_vec.resize(npar);
    force_vec.resize(npar, {});
    poteng_mat.resize(npar);
    force_cube.resize(npar);
    for (i=0; i<npar; i++) {
        poteng_mat[i].resize(npar);
        force_cube[i].resize(npar, {});
    }
}


/* -------------------------------------------------------------
   Compute the squared norm of a valarray 'array'
---------------------------------------------------------------- */

double ForceField::norm(std::valarray<double> array)
{
    double normsq = 0.;
    for (unsigned int i=0; i < array.size(); i++)
    {
        normsq += array[i] * array[i];
    }
    return normsq;
}


void ForceField::create_label_mapping()
{
    // find unique labels
    unique_labels = label1_vec;
    std::sort( unique_labels.begin(), unique_labels.end() );
    unique_labels.erase( std::unique( unique_labels.begin(),
                         unique_labels.end() ), unique_labels.end() );
    ntype = unique_labels.size();

    // map labels to types
    for (int i=0; i < ntype; i++){
        label2type[unique_labels[i]] = i;
    }
}
