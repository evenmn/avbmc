#include "fixed.h"
#include "../box.h"

Fixed::Fixed(class Box* box_in, const vec length_in)
    : Boundary(box_in)
{
    assert(length_in.n_cols == box->ndim);
    length = length_in;
    volume = prod(length);    
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
