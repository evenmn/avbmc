#include <iostream>
#include <cmath>
#include <string>
#include <valarray>
#include <memory>

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
    //boxes.push_back(box_in);
    dx = dx_in;
    label = "Trans   ";
}

/*
Trans::Trans(System* system_in, std::shared_ptr<Box> box_in, const double dx_in)
    : Moves(system_in)
{
    box = box_in;
    boxes.push_back(box_in);
    dx = dx_in;
    label = "Trans   ";
}
*/

/* ---------------------------------------------------
   Propose particle to move and new position randomly.
------------------------------------------------------ */

void Trans::perform_move()
{

    std::cout << "trans perform_move1" << std::endl;
    i = rng->next_int(box->npar); // particle to move
    std::cout << "trans perform_move2" << std::endl;
    double u0 = system->forcefield->comp_energy_par(box->particles, i);
    std::cout << "trans perform_move3" << std::endl;

    // move particle i
    std::valarray<double> dr(system->ndim);
    std::cout << "trans perform_move4" << std::endl;
    for(double &d : dr)
        d = 2 * (rng->next_double() - 0.5);
    std::cout << "trans perform_move5" << std::endl;
    pos_old = box->particles[i].r;
    std::cout << "trans perform_move6" << std::endl;
    box->particles[i].r += dx * dr;
    std::cout << "trans perform_move7" << std::endl;
    for (Particle 

    // compute new energy contribution from particle i
    double u1 = system->forcefield->comp_energy_par(box->particles, i);
    std::cout << "trans perform_move8" << std::endl;
    du = u1 - u0;
    std::cout << "trans perform_move9" << std::endl;
    box->poteng += du;
    std::cout << "trans perform_move10" << std::endl;
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
    box->particles[i].r = pos_old;
    box->poteng -= du;
}


/* ----------------------------------------------------------
   Update number of time this system size has occured if
   move was accepted
------------------------------------------------------------- */

void Trans::update_nsystemsize()
{
    if (box->npar + 1 > box->nsystemsize.size()) {
        box->nsystemsize.resize(box->npar + 1);
    }
    std::cout << "box->npar trans: " << box->npar << std::endl;
    std::cout << "box->nsystemsize.size(): " << box->nsystemsize.size() << std::endl;
    box->nsystemsize[box->npar] ++;
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
