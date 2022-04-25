#include <iostream>
#include <valarray>

#include "rungekutta4.h"
#include "../box.h"
#include "../particle.h"
#include "../forcefield/forcefield.h"


/* ------------------------
   Fourth order Runge-Kutta integration scheme
------------ */

RungeKutta4::RungeKutta4(Box* box_in, double dt_in)
    : Integrator(box_in, dt_in)
{
    dt2 = 0.5 * dt;
}


/* -----------------------
   Move particles to next step
-------------------- */

double RungeKutta4::next_step()
{
    unsigned int i;
    double energy;

    i = 0;
    energy = 0.;
    for (Particle &particle : box->particles) {
        // store old positions and velocities
        r_old = particle.r;
        v_old = particle.v;

        // calculate K1
        K1x = particle.v * dt2;
        K1v = particle.f * dt2 / box->mass[i];

        // calculate K2
        particle.r = r_old + K1x;
        particle.v = v_old + K1v;
        box->forcefield->comp_energy_par_force1(i, particle.f);
        K2x = particle.v * dt2;
        K2v = particle.f * dt2 / box->mass[i];

        // calculate K3
        particle.r = r_old + K2x;
        particle.v = v_old + K2v;
        box->forcefield->comp_energy_par_force1(i, particle.f);
        K3x = particle.v * dt;
        K3v = particle.f * dt / box->mass[i];

        // calculate K4
        particle.r = r_old + K3x;
        particle.v = v_old + K3v;
        box->forcefield->comp_energy_par_force1(i, particle.f);
        K4x = particle.v * dt;
        K4v = particle.f * dt / box->mass[i];
        
        // move
        particle.r = r_old + (K1x + 2 * (K2x + K3x) + K4x) / 6;
        particle.v = v_old + (K1v + 2 * (K2v + K3v) + K4v) / 6;
        energy += box->forcefield->comp_energy_par_force1(i, particle.f);
        i++;
    }
    return energy;
}
