#include "trans.h"
#include "../box.h"
#include "../particle.h"


Trans::Trans(class Box* box_in, double dx_in)
    : Moves(box_in) 
{
    dx = dx_in;
}

void Trans::perform_move()
{
    /* Propose naive translation move of particle i
     */

    // compute initial energy contribution from particle i
    i = box->rng->next_int(box->npar);
    double u0 = box->forcefield->comp_energy_par(box->particles, i);

    // move particle i
    std::valarray<double> dr(box->ndim);
    for(double &j : dr)
        j = rng->next_double();
    pos_old = box->particles[i]->r;
    box->particles[i]->r += 2 * (rng->next_double() - 0.5) * dx * dr;

    // compute new energy contribution from particle i
    double u1 = box->forcefield->comp_energy_par(box->particles, i);
    box->sampler->du = u1 - u0;
    box->poteng += box->sampler->du;
}

double Trans::accept(double /*temp*/, double /*chempot*/)
{
    /* Gives the ratio between the translation
     * probabilities Tij and Tji. For the 
     * naive case, Tij = Tji.
     */
    return 1.;
}

void Trans::reset()
{
    /* Set back to initial state if move is rejected
     */
    box->particles[i]->r = pos_old;
    box->poteng -= box->sampler->du;
}
