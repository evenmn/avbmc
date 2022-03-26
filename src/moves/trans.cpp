#include <iostream>
#include <cmath>
#include <string>
#include <valarray>

#include "trans.h"
#include "../box.h"
#include "../system.h"
#include "../particle.h"
#include "../rng/rng.h"
#include "../distance_manager.h"
#include "../boundary/boundary.h"
#include "../constraint/constraint.h"
#include "../forcefield/forcefield.h"


/*-----------------------------------------------------------------------------
  Constructor of naive translational move, with 'dx_in' being the maximum step
  length (in any direction)
------------------------------------------------------------------------------- */

Trans::Trans(System* system_in, Box* box_in, const double dx_in)
    : Moves(system_in)
{
    box = box_in;
    dx = dx_in;
    cum_time = 0.;
    label = "Trans   ";
}


/* ----------------------------------------------------------------------------
   Propose particle to move and new position randomly.
------------------------------------------------------------------------------- */

void Trans::perform_move()
{
    double u0, u1;

    i = rng->next_int(box->npar); // particle to move
    if (box->store_energy) {
        u0 = box->forcefield->poteng_vec[i];
    }
    else {
        u0 = box->forcefield->comp_energy_par_force0(i);
    }

    // move particle i
    std::valarray<double> dr(system->ndim);
    for(double &d : dr)
        d = 2 * (rng->next_double() - 0.5);
    pos_old = box->particles[i].r;
    box->particles[i].r += dx * dr;
    box->distance_manager->set();
    box->distance_manager->update_trans(i);

    // compute new energy contribution from particle i
    u1 = box->forcefield->comp_energy_par_force0(i);
    du = u1 - u0;
    box->poteng += du;
}


/* ----------------------------------------------------------------------------
   Compute acceptance probability of translational move
------------------------------------------------------------------------------- */

double Trans::accept(double temp, double /*chempot*/)
{
    bool constr_satis = true;
    for (Constraint* constraint : box->constraints) {
        constr_satis *= constraint->verify();
    }
    box->boundary->correct_position(i);
    return constr_satis * std::exp(-du/temp);
}


/* ----------------------------------------------------------------------------
   Set back to old state before move if move was rejected
------------------------------------------------------------------------------- */

void Trans::reset()
{
    box->distance_manager->reset();
    box->forcefield->reset();
    box->particles[i].r = pos_old;
    box->poteng -= du;
}


/* ----------------------------------------------------------------------------
   Update number of time this system size has occured if move was accepted
------------------------------------------------------------------------------- */

void Trans::update_nsystemsize()
{
    if (box->npar + 1 > box->nsystemsize.size()) {
        box->nsystemsize.resize(box->npar + 1);
    }
    box->nsystemsize[box->npar] ++;
}


/* ----------------------------------------------------------------------------
   Represent move in a clean way
------------------------------------------------------------------------------- */

std::string Trans::repr()
{
    std::string move_info;
    move_info += "Naive translational move\n";
    move_info += "    Maximum move length: " + std::to_string(dx) + "\n";
    return move_info;
}
