#include <iostream>
#include <string>
#include <valarray>
#include <cassert>

#include "periodic.h"
#include "../box.h"
#include "../system.h"
#include "../particle.h"


/* ----------------------
   Periodic boundary conditions class. Ensures that particles are always
   located inside the box. Potential energy and force calculations 
   always use the shortest distance between particles (might be through
   the wall). Velocities are not touched.
--------- */

Periodic::Periodic(Box* box_in, std::valarray<double> length_in)
    : Boundary(box_in)
{
    length = length_in;
    assert (length.size() == box->system->ndim);
    ndim = box->system->ndim;
    volume = 1.;
    for (double ld : length) {
        volume *= ld;
    }
    label = "Periodic";
}


/* -------------------
   Move particles to the other side of the box if they "cross" a wall.
   For this, we use the floor (modulus is an alternative)
------------ */

inline void Periodic::correct_position(unsigned int i)
{
    //box->particles[i].r = (box->particles[i].r + length) % length;
    for (unsigned int j=0; j<ndim; j++) {
        box->particles[i].r[j] -= floor(box->particles[i].r[j] / length[j]);
        //box->particles[i].r[j] = (box->particles[i].r[j] + length[j]) % length[j];
    }
}


/* -------------
   Find shortest distance between all particles. Often, this is through 
   the wall. For this, we use the round function on the relative
   coordinates.
------------------------- */

inline void Periodic::correct_distance(std::valarray<double> &dr)
{
    unsigned int i;
    for (i=0; i<ndim; i++) {
        dr[i] -= round(dr[i]/length[i]) * length[i];
    }
}
