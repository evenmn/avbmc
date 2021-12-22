#include <iostream>
#include <cmath>
#include <string>
#include <valarray>

#include "trans.h"
#include "../box.h"
#include "../system.h"
#include "../particle.h"
#include "../rng/rng.h"
#include "../forcefield/forcefield.h"


/*---------------------------------------------------
  Constructor of naive translational move, with
  'dx_in' being the maximum step length (in any
  direction)
----------------------------------------------------- */

Trans::Trans(System* system_in, Box* box_in, const double dx_in)
    : Moves(system_in)
{
    box = box_in;
    boxes.push_back(box_in);
    dx = dx_in;
    label = "Trans";
}


/* ---------------------------------------------------
   Propose particle to move and new position randomly.
------------------------------------------------------ */

void Trans::perform_move()
{
    i = rng->next_int(box->npar); // particle to move
    double u0 = system->forcefield->comp_energy_par(box->particles, i);

    // move particle i
    std::valarray<double> dr(system->ndim);
    for(double &d : dr)
        d = 2 * (rng->next_double() - 0.5);
    pos_old = box->particles[i]->r;
    box->particles[i]->r += dx * dr;

    // compute new energy contribution from particle i
    double u1 = system->forcefield->comp_energy_par(box->particles, i);
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
    std::cout << "reset1" << std::endl;
    box->particles[i]->r = pos_old;
    std::cout << "reset2" << std::endl;
    box->poteng -= du;
    std::cout << "reset3" << std::endl;
}


/* -----------------------------------------------------
   Represent move in a clean way
-------------------------------------------------------- */

std::string Trans::repr()
{
    std::string move_info;
    move_info += "Naive translational move\n";
    move_info += "    Maximum move length: " + std::to_string(dx) + "\n";
    return move_info;
}
