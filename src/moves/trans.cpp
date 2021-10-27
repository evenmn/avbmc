#include "trans.h"
#include "../box.h"


Trans::Trans(class Box* box_in, double dx_in)
    : Moves(box_in) 
{
    dx = dx_in;
}

void Trans::perform_move(const int i)
{
    /* Propose naive translation
     * move
     */

    // declare variables
    mat pos_copy = box->positions;
    rowvec dr(3), a0;
    double u0;

    // get acceleration and potential
    // energy of picked particle
    a0 = box->accelerations.row(i);
    u0 = box->potengs(i);

    // move particle i
    for(double &j : dr){
        j = rng->next_double();
    }
    dr *= dx / norm(dr);
    eps = 2 * (rng->next_double() - 0.5) * dr; 
    pos_copy.row(i) += eps;

    // compute new acceleration and potential energy
    u1 = box->forcefield->eval_acc_par(pos_copy, i, a1, true);
    box->sampler->da = a1 - a0;
    box->sampler->du = u1 - u0;
}

double Trans::accept()
{
    /* Gives the ratio between the translation
     * probabilities Tij and Tji. For the 
     * naive case, Tij = Tji.
     */
    return 1.;
}

void Trans::update_box(const int i)
{
    /* Update global (box) variables
     */
    box->positions.row(i) += eps;;
    box->accelerations.row(i) = a1;
    box->potengs(i) = u1;
    box->poteng += box->sampler->du;
}
