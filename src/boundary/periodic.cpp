/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.

  Author(s): Even M. Nordhagen
  Email(s): evenmn@mn.uio.no
  Date: 2022-06-03 (last changed 2022-06-03)
---------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------
   Periodic boundary conditions class. Ensures that particles are always
   located inside the box. Particles can go through the wall and end up on the
   opposite side of the box. Potential energy and force calculations always use
   the shortest distance between particles (might be through the wall),
   according to the minimum image convention. Velocities are left untouched.
---------------------------------------------------------------------------- */

#include <iostream>
#include <string>
#include <valarray>
#include <cassert>
#include <numeric>  // accumulate
#include <cmath>    // fmod

#include "periodic.h"
#include "../box.h"
#include "../system.h"
#include "../particle.h"


/* ----------------------------------------------------------------------------
   Periodic boundary conditions class constructor.
---------------------------------------------------------------------------- */

Periodic::Periodic(Box* box_in, std::valarray<double> length_in)
    : Boundary(box_in), length(length_in)
{
    assert (length.size() == box->system->ndim);
    ndim = box->system->ndim;
    //volume = std::accumulate(length.begin(), length.end(), 1, std::multiplies<double>);
    volume = 1.;
    for (double ld : length) {
        volume *= ld;
    }
    label = "Periodic";
}


/* ----------------------------------------------------------------------------
   Move particles to the other side of the box if they "cross" a wall.
   For this, we use the floor (modulus is an alternative)
---------------------------------------------------------------------------- */

inline void Periodic::correct_position(unsigned int i)
{
    unsigned int j;

    for (j=0; j<ndim; j++) {
        box->particles[i].r[j] = std::fmod(box->particles[i].r[j], length[j]);
    }
}


/* ----------------------------------------------------------------------------
   Find shortest distance between all particles. Often, this is through 
   the wall. For this, we use the round function on the relative
   coordinates.
---------------------------------------------------------------------------- */

inline void Periodic::correct_distance(std::valarray<double> &dr)
{
    unsigned int i;

    for (i=0; i<ndim; i++) {
        dr[i] -= round(dr[i]/length[i]) * length[i];
    }
}
