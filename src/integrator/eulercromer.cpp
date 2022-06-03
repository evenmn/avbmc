/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.

  Author(s): Even M. Nordhagen
  Email(s): evenmn@mn.uio.no
  Date: 2022-06-03 (last changed 2022-06-03)
---------------------------------------------------------------------------- */

#include <iostream>
#include <valarray>

#include "eulercromer.h"
#include "../box.h"
#include "../particle.h"
#include "../forcefield/forcefield.h"


/* ----------------------------------------------------------------------------
   Euler Cromer integration class. Is proven to be stable for
   oscillatory problems in mechanics, see A. Cromer, American Journal of Physics 49, 455 (1981)
------------------------------------------------------------------------------- */

EulerCromer::EulerCromer(Box* box_in, double dt_in)
    : Integrator(box_in, dt_in)
{}


/* ----------------------------------------------------------------------------
   Move to next step
------------------------------------------------------------------------------- */

double EulerCromer::next_step()
{
    unsigned int i;
    double energy;

    i = 0;
    energy = 0.;
    for (Particle &particle : box->particles) {
        energy += box->forcefield->comp_energy_par_force1(i, particle.f);
        particle.v += particle.f * dt / particle.mass;
        particle.r += particle.v * dt;
        i++;
    }
    return energy;
}
