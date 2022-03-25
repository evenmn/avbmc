#include <iostream>
#include <valarray>

#include "euler.h"
#include "../box.h"
#include "../particle.h"
#include "../forcefield/forcefield.h"


/* ------------------
   Euler integration class
---------------------- */

Euler::Euler(Box* box_in, double dt_in)
    : Integrator(box_in, dt_in)
{}


/* -----------------------------
   Move to next step
------------------ */

double Euler::next_step()
{
    unsigned int i;
    double energy;

    i = 0;
    energy = 0.;
    for (Particle &particle : box->particles) {
        energy + = box->forcefield->comp_energy_par_force1(i, particle.f);
        particle.r += particle.v * dt;
        particle.v += particle.f * dt / box->mass[i];
        i++;
    }
    return energy;
}
