#include <iostream>
#include <string>
#include <vector>
#include <valarray>

#include "forcefield.h"
#include "../system.h"
#include "../particle.h"

ForceField::ForceField(System* system_in)
{
    system = system_in;
    ntype = nline = 0;
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
