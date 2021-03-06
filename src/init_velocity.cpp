/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.

  Author(s): Even M. Nordhagen
  Email(s): evenmn@mn.uio.no
  Date: 2022-06-03 (last changed 2022-06-03)
---------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------
  Get initial particle velocities. NB: This interface is outdated and should
  be updated when molecular dynamics is implemented.
---------------------------------------------------------------------------- */

#include <iostream>
#include <valarray>

#include "init_velocity.h"
#include "rng/rng.h"


/* ----------------------------------------------------------------------------
   Zero velocity constructor.
---------------------------------------------------------------------------- */

Zero::Zero() : Velocity() {}

std::vector<std::valarray<double> > Zero::get_velocity(unsigned int npar,
    unsigned int ndim)
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


/* ----------------------------------------------------------------------------
   Random normal velocity constructor.
---------------------------------------------------------------------------- */

std::vector<std::valarray<double> > Gauss::get_velocity(unsigned int npar,
    unsigned int ndim)
{
    unsigned int i, j;
    std::vector<std::valarray<double> >  velocities;
    for (i=0; i<npar; i++) {
        std::valarray<double> velocity(0., ndim);
        velocities.push_back(velocity);
        for (j=0; j<ndim; j++) {
            velocities[i][j] = rng->next_gaussian(mean, var);
        }
    }
    return velocities;
}


/* ----------------------------------------------------------------------------
   Initialize system at some temperature by assigning random velocities drawn
   from a normal distribution. The variance of the normal distribution is 
   chosen according to the equipartition theorem.
---------------------------------------------------------------------------- */

Temp::Temp(RandomNumberGenerator* rng_in, double temp_in)
    : Gauss(rng_in, 0., sqrt(temp_in)) {}
