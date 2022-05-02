#include <iostream>
#include <valarray>

#include "velocityverlet.h"
#include "../box.h"
#include "../particle.h"
#include "../forcefield/forcefield.h"


/* ----------------------------------------------------------------------------
   Velocity Verlet integrator class
------------------------------------------------------------------------------- */

VelocityVerlet::VelocityVerlet(Box* box_in, double dt_in)
    : Integrator(box_in, dt_in)
{
    dt2 = 0.5 * dt;
    ddt2 = dt * dt2;
}


/* ----------------------------------------------------------------------------
   Move to next step
------------------------------------------------------------------------------- */

double VelocityVerlet::next_step()
{
    unsigned int i;
    double energy;
    std::valarray<double> f_old;

    i = 0;
    energy = 0.;
    for (Particle &particle : box->particles) {
        f_old = particle.f;
        particle.r += particle.v * dt + particle.f * ddt2 / particle.mass;
        energy += box->forcefield->comp_energy_par_force1(i, particle.f);
        particle.v += (particle.f + f_old) * dt2 / particle.mass;
        i++;
    }
    return energy;
}
