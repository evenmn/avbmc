#include <iostream>
#include <cmath>
#include <valarray>

#include "trans.h"
#include "../box.h"
#include "../particle.h"


/*---------------------------------------------------
  Constructor of naive translational move, with
  'dx_in' being the maximum step length (in any
  direction)
----------------------------------------------------- */

Trans::Trans(Box* box_in, double dx_in)
    : Moves(box_in) 
{
    dx = dx_in;
}


/* ---------------------------------------------------
   Propose particle to move and new position randomly.
------------------------------------------------------ */

void Trans::perform_move()
{
    i = box->rng->next_int(box->npar); // particle to move
    double u0 = box->forcefield->comp_energy_par(box->particles, i);

    // move particle i
    std::valarray<double> dr(box->ndim);
    for(double &j : dr)
        j = rng->next_double();
    pos_old = box->particles[i]->r;
    box->particles[i]->r += 2 * (rng->next_double() - 0.5) * dx * dr;

    // compute new energy contribution from particle i
    double u1 = box->forcefield->comp_energy_par(box->particles, i);
    du = u1 - u0;
    box->poteng += du;
}


/* -----------------------------------------------------
   Compute acceptance probability of translational 
   move
-------------------------------------------------------- */

double Trans::accept(double temp, double /*chempot*/)
{
    return std::exp(-du/temp);
}


/* -----------------------------------------------------
   Set back to old state before move if move was
   rejected
-------------------------------------------------------- */

void Trans::reset()
{
    box->particles[i]->r = pos_old;
    box->poteng -= du;
}
