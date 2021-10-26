#include "transmh.h"
#include "../box.h"


TransMH::TransMH(class Box* box_in, double dx_in, double Ddt_in)
    : Moves(box_in) 
{
    dx = dx_in;
    Ddt = Ddt_in;
}

mat TransMH::perform_move(const int i)
{
    /* Propose translation move
     * based on Metropolis-Hastings
     * criterion
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
        j = rng->next_gaussian();
    }
    dr *= dx / norm(dr);
    eps = Ddt * a0 + dr;
    pos_copy.row(i) += eps;

    // check Stillinger criterion with new position
    bool stillinger_accepted = check_stillinger();
    if(stillinger_accepted){
        // compute new acceleration and potential energy
        u1 = box->forcefield->eval_acc_par(pos_copy, i, a1, true);
        box->sampler->da = a1 - a0;
        box->sampler->du = u1 - u0;
        return pos_copy;
    }
    else{
        box->sampler->da = zeros(3);
        box->sampler->du = 0.;
        return box->positions;
    }
}

double TransMH::accept()
{
    /* Gives the ratio between the translation
     * probabilities Tij and Tji. This is the ratio
     * between two Green's functions. 1 comes from
     * when da=0
     */
    return exp(0.5 * as_scalar(box->sampler->da*eps)) + 1;
}

void TransMH::update_box(const int i)
{
    /* If the move is accepted, we want to update
     * the position matrix, acceleration matrix
     * and potential energy vector/double since
     * this is already computed.
     */
    box->positions.row(i) += eps;;
    box->accelerations.row(i) = a1;
    box->potengs(i) = u1;
    box->poteng += box->sampler->du;
}
