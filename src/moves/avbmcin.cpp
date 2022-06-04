/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.

  Author(s): Even M. Nordhagen
  Email(s): evenmn@mn.uio.no
  Date: 2022-06-03 (last changed 2022-06-03)
---------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------
   Aggregation-volume-biased Monte Carlo (AVBMC) insertion move. This is a
   fictitious inter-box move, and works in the grand-canonical essemble only.
   To maintain detailed balance, this move always have to be used in
   conjunction with an AVBMC deletion move, with the same radius of the 
   bonded region and the same move probability. This is ensured when applying
   the 'AVBMC' move. The AVBMC moves were first proposed by Chen (2000).
---------------------------------------------------------------------------- */

#include <iostream>
#include <cmath>
#include <valarray>
#include <vector>

#include "avbmcin.h"
#include "../box.h"
#include "../system.h"
#include "../particle.h"
#include "../rng/rng.h"
#include "../sampler/sampler.h"
#include "../boundary/boundary.h"
#include "../forcefield/forcefield.h"
#include "../distance_manager.h"


/* ----------------------------------------------------------------------------
   AVBMCIn constructor
---------------------------------------------------------------------------- */

AVBMCIn::AVBMCIn(System* system_in, Box* box_in, const std::string &particle_in,
    const double r_below_in, const double r_above_in, bool energy_bias_in)
    : Moves(system_in)
{
    box = box_in;
    energy_bias = energy_bias_in;
    r_below = r_below_in;
    r_above = r_above_in;
    r_abovesq = r_above * r_above;
    r_belowsq = r_below * r_below;
    v_in = 1.;

    particle_label = particle_in;
    particle_type = box->forcefield->label2type.at(particle_in);
    label = "AVBMCIn ";
}


/* ----------------------------------------------------------------------------
   Insert molecule into the bonded region 
---------------------------------------------------------------------------- */

void AVBMCIn::perform_move()
{
    unsigned int i;
    std::valarray<double> dr(system->ndim);

    // pick target particle and count number of neighbors
    i = box->typeidx[particle_type][rng->next_int(box->npartype[particle_type])];
    std::vector<unsigned int> neigh_listi = box->build_neigh_list(i, r_abovesq);
    n_in = neigh_listi.size();

    // construct new particle close to target particle
    dr = insertion_position(false);
    box->add_particle(particle_label, box->particles[i].r + dr);

    // compute du (also bake this into add_particle?)
    du = box->forcefield->comp_energy_par_force0(box->npar - 1);
    box->poteng += du;
}


/* ----------------------------------------------------------------------------
   Get acceptance probability of move
---------------------------------------------------------------------------- */

double AVBMCIn::accept(double temp, double chempot)
{
    double dw = system->sampler->w(box->npar) - system->sampler->w(box->npar-1);
    double prefactor = v_in * box->npar / ((n_in + 1) * (box->npar + 1));
    return prefactor * std::exp(-(du-chempot+dw)/temp);
}


/* ----------------------------------------------------------------------------
   Set back to old state before move is rejected
---------------------------------------------------------------------------- */

void AVBMCIn::reset()
{
    box->rm_particle(box->npar-1);
    box->poteng -= du;
}


/* ----------------------------------------------------------------------------
   Update number of time this system size has occured if move was accepted
---------------------------------------------------------------------------- */

void AVBMCIn::update_size_histogram()
{
    box->update_size_histogram();
}


/* ----------------------------------------------------------------------------
   Represent move in a clean way
---------------------------------------------------------------------------- */

std::string AVBMCIn::repr()
{
    std::string move_info;
    move_info += "AVBMC particle insertion move\n";
    move_info += "    Radius of outer sphere: " + std::to_string(r_above) + "\n";
    move_info += "    Radius of inner sphere: " + std::to_string(r_below) + "\n";
    move_info += "    Label of inserted atom: " + particle_label + "\n";
    return move_info;
}
