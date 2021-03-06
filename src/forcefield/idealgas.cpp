/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.

  Author(s): Even M. Nordhagen
  Email(s): evenmn@mn.uio.no
  Date: 2022-06-03 (last changed 2022-06-03)
---------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------
  The ideal gas forcefield can be useful when sampling a gas reservoir, for
  instance. No force will act between particles, and all moves will be 
  accepted. It is necessary to specify the expected elements in order for
  the forcefield class to create type mapping like for other forcefields.
---------------------------------------------------------------------------- */

#include <iostream>
#include <valarray>

#include "idealgas.h"
#include "../box.h"
#include "../system.h"


/* ----------------------------------------------------------------------------
   This is the default constructor when a parameter file 'params' is given.
---------------------------------------------------------------------------- */

IdealGas::IdealGas(Box* box_in, const std::vector<std::string> &elements_in)
    : ForceField(box_in)
{
    label = "Ideal gas";
    label1_vec = elements_in;
    create_label_mapping();
    init_ntype();
}



/* ----------------------------------------------------------------------------
   Compute energy contribution from a particle 'i' without using neighbor
   lists.
---------------------------------------------------------------------------- */

double IdealGas::comp_energy_par_neigh0_eng0(const unsigned int /*i*/,
    std::valarray<double> &force, const bool /*comp_force*/)
{
    force.resize(box->system->ndim, 0.);
    return 0.;
}


/* ----------------------------------------------------------------------------
   Compute energy contribution from a particle 'i'
---------------------------------------------------------------------------- */

double IdealGas::comp_energy_par_neigh1_eng0(const unsigned int /*i*/,
    std::valarray<double> &force, const bool /*comp_force*/)
{
    force.resize(box->system->ndim, 0.);
    return 0.;
}


/* ----------------------------------------------------------------------------
   Compute energy contribution from a particle 'i'
---------------------------------------------------------------------------- */

double IdealGas::comp_energy_par_neigh1_eng1(const unsigned int /*i*/,
    std::valarray<double> &force, const bool /*comp_force*/)
{
    force.resize(box->system->ndim, 0.);
    return 0.;
}
