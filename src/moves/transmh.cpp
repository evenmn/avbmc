#include <iostream>
#include <cmath>
#include <string>
#include <valarray>

#include "transmh.h"
#include "../box.h"
#include "../system.h"
#include "../particle.h"
#include "../rng/rng.h"
#include "../boundary/boundary.h"
#include "../forcefield/forcefield.h"


/*---------------------------------------------------
  Constructor of biased translational move, with
  'dx_in' being the maximum step length (in any
  direction) and 'Ddt_in' being the biasing 
  parameter.
----------------------------------------------------- */

TransMH::TransMH(System* system_in, Box* box_in, const double dx_in,
                 const double Ddt_in)
    : Moves(system_in)
{
    box = box_in;
    dx = dx_in;
    Ddt = Ddt_in;
    label = "TransMH ";
}


/* ---------------------------------------------------
   Propose particle to move and new position randomly.
------------------------------------------------------ */

void TransMH::perform_move()
{
    i = rng->next_int(box->npar); // particle to move
    std::valarray<double> f0;
    double u0 = system->forcefield->comp_energy_par(box->particles, i, f0, true);

    // move particle i
    std::valarray<double> dr(system->ndim);
    for(double &d : dr)
        d = 2 * (rng->next_double() - 0.5);
    dr *= dx / norm(dr);
    eps = Ddt * f0 + dr;
    box->particles[i].r += eps;

    // compute new energy contribution from particle i
    std::valarray<double> f1;
    double u1 = system->forcefield->comp_energy_par(box->particles, i, f1, true);
    du = u1 - u0;
    df = f1 - f0;
    box->poteng += du;
}


/* -----------------------------------------------------
   Compute acceptance probability of translational
   move
-------------------------------------------------------- */

double TransMH::accept(double temp, double /*chempot*/)
{
    double dot = 0.;
    for (unsigned int j = 0; j < system->ndim; j++) {
        dot += df[i] * eps[i];
    }
    return (std::exp(0.5 * dot) + 1) * std::exp(-du/temp) * box->boundary->correct_position();
}


/* -----------------------------------------------------
   Set back to old state before move if move was
   rejected
-------------------------------------------------------- */

void TransMH::reset()
{
    box->particles[i].r -= eps;
    box->poteng -= du;
}


/* ----------------------------------------------------------
   Update number of time this system size has occured if
   move was accepted
------------------------------------------------------------- */

void TransMH::update_nsystemsize()
{
    if (box->npar + 1 > box->nsystemsize.size()) {
        box->nsystemsize.resize(box->npar + 1);
    }
    box->nsystemsize[box->npar] ++;
}


/* -----------------------------------------------------
   Represent move in a clean way
-------------------------------------------------------- */

std::string TransMH::repr()
{
    std::string move_info;
    move_info += "Biased translational move\n";
    move_info += "    Maximum move length: " + std::to_string(dx) + "\n";
    move_info += "    Biasing parameter:   " + std::to_string(Ddt) + "\n";
    return move_info;
}
