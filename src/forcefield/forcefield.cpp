#include <iostream>
#include <string>
#include <vector>
#include <valarray>

#include "forcefield.h"
#include "../box.h"
#include "../system.h"
#include "../particle.h"


/* ----------------------------------------------------------------------------
   Forcefield base class. Making 'comp_energy_par' point to the correct energy
   calculation function
------------------------------------------------------------------------------- */

ForceField::ForceField(Box* box_in)
{
    box = box_in;
    ntype = nline = 0;

    // let comp_energy_par point to the desired energy function
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


/* ----------------------------------------------------------------------------
   Compute energy contribution from a particle 'i', when force is not needed.
------------------------------------------------------------------------------- */

double ForceField::comp_energy_par_force0(const unsigned int i)
{
    std::valarray<double> force;
    return (this->*(this->comp_energy_par))(i, force, false);
}


/* ----------------------------------------------------------------------------
   Compute energy contribution from a particle 'i', when force is needed.
------------------------------------------------------------------------------- */

double ForceField::comp_energy_par_force1(const unsigned int i,
                                          std::valarray<double> &force)
{
    return (this->*(this->comp_energy_par))(i, force, true);
}


/* ----------------------------------------------------------------------------
   Initialize energy and force matrices
------------------------------------------------------------------------------- */

void ForceField::initialize()
{
    unsigned int npar, i;

    // initialize matrices
    npar = box->npar;
    poteng_vec.resize(npar);
    force_vec.resize(npar, {});
    poteng_mat.resize(npar);
    force_cube.resize(npar);
    for (i=0; i<npar; i++) {
        poteng_mat[i].resize(npar);
        force_cube[i].resize(npar, {});
    }
    poteng_vec_old = poteng_vec;
    force_vec_old = force_vec;
    poteng_mat_old = poteng_mat;
    force_cube_old = force_cube;

    // fill matrices
    //comp_energy_all(
    //for (i=0; i<npar; i++) {

}


/* ----------------------------------------------------------------------------
   Store current energy and force matrices
------------------------------------------------------------------------------- */

void ForceField::set()
{
    poteng_vec_old = poteng_vec;
    force_vec_old = force_vec;
    poteng_mat_old = poteng_mat;
    force_cube_old = force_cube;
}


/* ----------------------------------------------------------------------------
   Reset energy and force matrices
------------------------------------------------------------------------------- */

void ForceField::reset()
{
    poteng_vec = poteng_vec_old;
    force_vec = force_vec_old;
    poteng_mat = poteng_mat_old;
    force_cube = force_cube_old;
} 


/* ----------------------------------------------------------------------------
   Compute the squared norm of a valarray 'array'
------------------------------------------------------------------------------- */

double ForceField::norm(std::valarray<double> array)
{
    unsigned int i;
    double normsq;
    
    normsq = 0.;
    for (i=0; i < array.size(); i++)
    {
        normsq += array[i] * array[i];
    }
    return normsq;
}


/* ----------------------------------------------------------------------------
   Create mapping from label to type
------------------------------------------------------------------------------- */

void ForceField::create_label_mapping()
{
    unsigned int i;

    unique_labels = label1_vec;
    std::sort( unique_labels.begin(), unique_labels.end() );
    unique_labels.erase( std::unique( unique_labels.begin(),
                         unique_labels.end() ), unique_labels.end() );
    ntype = unique_labels.size();

    // map labels to types
    for (i=0; i < ntype; i++){
        label2type[unique_labels[i]] = i;
    }
}
