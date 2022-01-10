#include <iostream>
#include <vector>
#include <valarray>

#include "forcefield.h"
#include "../system.h"
#include "../particle.h"

ForceField::ForceField(System* system_in)
{
    system = system_in;
}


/* -------------------------------------------------------
   Compute energy contribution from a particle 'i'
---------------------------------------------------------- */

double ForceField::comp_energy_par(std::vector<Particle> particles, const int i)
{
    std::valarray<double> force;
    return (comp_energy_par(particles, i, force, false));
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
