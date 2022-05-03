#include <iostream>
#include <valarray>

#include "init_velocity.h"
#include "rng/rng.h"


/* ----------------------------------------------------------------------------
   Zero velocity constructor.
------------------------------------------------------------------------------- */

Zero::Zero() : Velocity() {}

std::vector<std::valarray<double> > Zero::get_velocity(unsigned int npar, unsigned int ndim)
{
    unsigned int i;

    std::vector<std::valarray<double> > velocities;
    for (i=0; i<npar; i++) {
        std::valarray<double> velocity(0., ndim);
        velocities.push_back(velocity);
    }
    return velocities;
}

Gauss::Gauss(RandomNumberGenerator* rng_in, double mean_in, double var_in) 
    : Velocity()
{
    rng = rng_in;
    mean = mean_in;
    var = var_in;
}

std::vector<std::valarray<double> > Gauss::get_velocity(unsigned int npar, unsigned int ndim)
{
    unsigned int i, j;
    std::vector<std::valrray<double> >  velocities;
    for (i=0; i<npar; i++) {
        std::valarray<double> velocity(0., ndim);
        velocities.push_back(velocity);
        for (j=0; j<ndim; j++) {
            velocities[i][j] = rng->next_gaussian(mean, var);
        }
    }
    return velocities;
}

Temp::Temp(RandomNumberGenerator* rng_in, double temp_in)
    : Gauss(rng_in, 0., sqrt(temp_in)) {}
