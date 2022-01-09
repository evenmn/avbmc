#include <iostream>
#include <valarray>

#include "../box.h"
#include "../particle.h"
#include "fixed.h"


/* ----------------------------------------------------------------------------
   Fixed boundary constructor
------------------------------------------------------------------------------- */

Fixed::Fixed(Box* box_in, std::valarray<double> length_in)
    : Boundary(box_in)
{
    length = length_in;
    boxdim = length.size();
    volume = 1.;
    for (double ld : length) {
        volume *= ld;
    }
}


/* ----------------------------------------------------------------------------
   Check if particles are outside of box. If this is the case, 'false' is
   returned and the Monte Carlo move is rejected. Molecular dynamics
   simulations should crash if that's the case.
------------------------------------------------------------------------------- */

bool Fixed::correct_position()
{
    for (Particle particle : box->particles) {
        for (unsigned int i = 0; i < boxdim; i++) {
            if (particle.r[i] < 0. || particle.r[i] > length[i]) {
                return false;
            }
        }
    }
    return true;
}


/* ----------------------------------------------------------------------------
   Return the volume of the box. This is simply the product of the lengths,
   which is already calculated in the constructor.
------------------------------------------------------------------------------- */

double Fixed::comp_volume()
{
    return volume;
}
