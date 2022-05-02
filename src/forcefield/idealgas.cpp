#include <iostream>
#include <valarray>

#include "idealgas.h"
#include "../box.h"
#include "../system.h"


/* ----------------------------------------------------------------------------
   This is the default constructor when a parameter file 'params' is given.
------------------------------------------------------------------------------- */

IdealGas::IdealGas(Box* box_in, const std::vector<std::string> &elements_in)
    : ForceField(box_in)
{
    label = "Ideal gas";
    label1_vec = elements_in;
    create_label_mapping();
}



/* ----------------------------------------------------------------------------
   Compute energy contribution from a particle 'i' without using neighbor
   lists.
------------------------------------------------------------------------------- */

double IdealGas::comp_energy_par_neigh0_eng0(const unsigned int /*i*/,
    std::valarray<double> &force, const bool /*comp_force*/)
{
    force.resize(box->system->ndim, 0.);
    return 0.;
}


/* ----------------------------------------------------------------------------
   Compute energy contribution from a particle 'i'
------------------------------------------------------------------------------- */

double IdealGas::comp_energy_par_neigh1_eng0(const unsigned int /*i*/,
    std::valarray<double> &force, const bool /*comp_force*/)
{
    force.resize(box->system->ndim, 0.);
    return 0.;
}


/* ----------------------------------------------------------------------------
   Compute energy contribution from a particle 'i'
------------------------------------------------------------------------------- */

double IdealGas::comp_energy_par_neigh1_eng1(const unsigned int /*i*/,
    std::valarray<double> &force, const bool /*comp_force*/)
{
    force.resize(box->system->ndim, 0.);
    return 0.;
}
