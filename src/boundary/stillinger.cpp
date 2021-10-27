#include "stillinger.h"
#include "../box.h"

Stillinger::Stillinger(class Box* box_in, double r_c_in)
    : Boundary(box_in)
{
    r_c = r_c_in;
    v_c = 4 * datum::pi * pow(r_c, 3) / 3;
}

bool Stillinger::correct_position(mat &pos)
{
    /* Check if Stillinger criterion is satisfied.
     * If it is not satisfied, Monte Carlo moves
     * should be rejected. Molecular dynamics
     * simulations should abort in the same case.
     */

    // check if criterion is satisfied
    return true;
}

bool Stillinger::correct_velocity(mat &vel)
{
    /* For Stillinger cluster criterion, the
     * velocity does not need to be corrected
     */
    return true;
}

bool Stillinger::correct_distance(mat &dist)
{
    /* For Stillinger cluster criterion, the
     * distance does not need to be corrected
     */
    return true;
}

double Stillinger::comp_volume()
{
    /* Get volume of cluster, but ignoring
     * the factor pi/3. This is done by first
     * overestimating the volume as N * v_stillinger,
     * and then subtracting the volume of overlapping
     * sphere caps.
     */

    double v_overest = box->npar * v_c;
    
    mat distance_mat_sqrd = box->forcefield->distance_mat;

    // compute height of all caps. This might be time-consuming
    mat h = (real(sqrtmat(distance_mat_sqrd)) / 2. - r_c);
    double v_caps = 2 * sum(sum(powmat(h, 2)%(-h + 3 * r_c)));

    return v_overest - v_caps;
}
