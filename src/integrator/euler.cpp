/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.

  Author(s): Even M. Nordhagen
  Email(s): evenmn@mn.uio.no
  Date: 2022-06-03 (last changed 2022-06-03)
---------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------
  Forward Euler integrator. Known to be unstable for most applications and
  should only be used for teaching purposes.
---------------------------------------------------------------------------- */

#include <iostream>
#include <valarray>

#include "euler.h"
#include "../box.h"
#include "../particle.h"
#include "../forcefield/forcefield.h"


/* ----------------------------------------------------------------------------
   Euler integration class
---------------------------------------------------------------------------- */

Euler::Euler(Box* box_in, double dt_in)
    : Integrator(box_in, dt_in)
{}


/* ----------------------------------------------------------------------------
   Move to next step
---------------------------------------------------------------------------- */

double Euler::next_step()
{
    unsigned int i;
    double energy;

    i = 0;
    energy = 0.;
    for (Particle &particle : box->particles) {
        energy += box->forcefield->comp_energy_par_force1(i, particle.f);
        particle.r += particle.v * dt;
        particle.v += particle.f * dt / particle.mass;
        i++;
    }
    return energy;
}
