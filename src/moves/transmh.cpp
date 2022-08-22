/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.

  Author(s): Even M. Nordhagen
  Email(s): evenmn@mn.uio.no
  Date: 2022-06-03 (last changed 2022-06-03)
---------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------
  Biased translational move. Random particle is moved a random distance biased
  in the direction of the force. It is accepted according to the Boltzmann
  distribution. 
---------------------------------------------------------------------------- */

#include <iostream>
#include <cmath>
#include <string>
#include <valarray>

#include "transmh.h"
#include "../box.h"
#include "../system.h"
#include "../particle.h"
#include "../distance_manager.h"
#include "../rng/rng.h"
#include "../boundary/boundary.h"
#include "../forcefield/forcefield.h"


/*-----------------------------------------------------------------------------
  Constructor of biased translational move, with 'dx_in' being the maximum step
  length (in any direction) and 'Ddt_in' being the biasing parameter.
---------------------------------------------------------------------------- */

TransMH::TransMH(System* system_in, Box* box_in, const double dx_in,
                 const double Ddt_in, const std::string &element_in)
    : Moves(system_in)
{
    box = box_in;
    dx = dx_in;
    Ddt = Ddt_in;
    cum_time = 0.;
    naccept = ndrawn = 0;
    if (element_in.empty()) {
        element_spec = false;
    }
    else {
        element_type = box->forcefield->label2type.at(element_in);
        element_spec = true;
    }
    label = "TransMH ";
}


/* ----------------------------------------------------------------------------
   Propose particle to move and new position randomly.
---------------------------------------------------------------------------- */

void TransMH::perform_move()
{
    double u0, u1;
    std::valarray<double> f0, f1;

    if (element_spec) {
        i = box->typeidx[element_type][rng->next_int(box->npartype[element_type])];
    }
    else {
        i = rng->next_int(box->npar);
    }

    if (box->store_energy) {
        u0 = box->forcefield->poteng_vec[i];
        f0 = box->forcefield->force_vec[i];
    }
    else {
        u0 = box->forcefield->comp_energy_par_force1(i, f0);
    }

    // move particle i
    std::valarray<double> dr(system->ndim);
    for(double &d : dr)
        d = 2 * (rng->next_double() - 0.5);
    dr *= dx / norm(dr);
    eps = Ddt * f0 + dr;
    box->particles[i].r += eps;
    box->distance_manager->set();
    box->distance_manager->update_trans(i);

    // compute new energy contribution from particle i
    u1 = box->forcefield->comp_energy_par_force1(i, f1);
    du = u1 - u0;
    df = f1 - f0;
    box->poteng += du;
}


/* ----------------------------------------------------------------------------
   Compute acceptance probability of translational move
---------------------------------------------------------------------------- */

double TransMH::accept(double temp, double /*chempot*/)
{
    double dot;
    unsigned int j;
    dot = 0.;
    for (j = 0; j < system->ndim; j++) {
        dot += df[j] * eps[j];
    }
    box->boundary->correct_position(i);
    return (std::exp(0.5 * dot) + 1) * std::exp(-du/temp);
}


/* ----------------------------------------------------------------------------
   Set back to old state before move if move was rejected
---------------------------------------------------------------------------- */

void TransMH::reset()
{
    box->distance_manager->reset();
    box->forcefield->reset();
    box->particles[i].r -= eps;
    box->poteng -= du;
}


/* ----------------------------------------------------------------------------
   Update number of time this system size has occured if move was accepted
---------------------------------------------------------------------------- */

void TransMH::update_size_histogram()
{
    box->update_size_histogram();
}


/* ----------------------------------------------------------------------------
   Represent move in a clean way
---------------------------------------------------------------------------- */

std::string TransMH::repr()
{
    std::string move_info;
    move_info += "Biased translational move\n";
    move_info += "    Maximum move length: " + std::to_string(dx) + "\n";
    move_info += "    Biasing parameter:   " + std::to_string(Ddt) + "\n";
    return move_info;
}
