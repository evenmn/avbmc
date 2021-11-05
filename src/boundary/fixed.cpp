#include "fixed.h"
#include "../box.h"

Fixed::Fixed(class Box* box_in, const double lx_in)
    : Boundary(box_in)
{
    assert(box->ndim == 1);
    lx = lx_in;
    volume = lx;    
}

Fixed::Fixed(class Box* box_in, const double lx_in, const double ly_in)
    : Boundary(box_in)
{
    assert(box->ndim == 2);
    lx = lx_in;
    ly = ly_in;
    volume = lx * ly;    
}

Fixed::Fixed(class Box* box_in, const double lx_in, const double ly_in, const double lz_in)
    : Boundary(box_in)
{
    assert(box->ndim == 3);
    lx = lx_in;
    ly = ly_in;
    lz = lz_in;
    volume = lx * ly * lz;    
}

bool Fixed::correct_position(mat &pos)
{
    /* Check if particles are outside box.
     * If they are, Monte Carlo moves
     * should be rejected. Molecular dynamics
     * simulations should abort in the same case.
     */

    // check if criterion is satisfied
    return true;
}

bool Fixed::correct_velocity(mat &vel)
{
    /* For fixed boundaries, the
     * velocity does not need to be corrected
     */
    return true;
}

bool Fixed::correct_distance(mat &dist)
{
    /* For fixed boundaries, the
     * distance does not need to be corrected
     */
    return true;
}

double Fixed::comp_volume()
{
    /* Get volume of box
     */
    return volume;
}
