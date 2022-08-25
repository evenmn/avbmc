/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.

  Author(s): Even M. Nordhagen
  Email(s): evenmn@mn.uio.no
  Date: 2022-06-03 (last changed 2022-08-19)
---------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------
   Aggregation-volume-biased Monte Carlo (AVBMC) deletion move. This is a
   fictitious inter-box move, and works in the grand-canonical essemble only.
   To maintain detailed balance, this move always have to be used in
   conjunction with an AVBMC insertion move, with the same radius of the 
   bonded region and the same move probability. This is ensured when applying
   the 'AVBMC' move. The AVBMC moves were first proposed by Chen (2000).
---------------------------------------------------------------------------- */

#include <iostream>
#include <cmath>
#include <vector>
#include <random>

#include "avbmcout.h"
#include "../box.h"
#include "../system.h"
#include "../particle.h"
#include "../rng/rng.h"
#include "../sampler/sampler.h"
#include "../boundary/boundary.h"
#include "../forcefield/forcefield.h"
#include "../distance_manager.h"


/* ----------------------------------------------------------------------------
   AVBMCOut constructor
---------------------------------------------------------------------------- */

AVBMCOut::AVBMCOut(System* system_in, Box* box_in,
    const std::string &particle_in, const double r_above_in, bool energy_bias_in) 
    : Moves(system_in), particle_label(particle_in)
{
    box = box_in;
    energy_bias = energy_bias_in;
    r_above = r_above_in;
    r_abovesq = r_above * r_above;
    v_in = 1.; // 4 * pi * std::pow(r_above, 3)/3;
    cum_time = 0.;
    n_in = naccept = ndrawn = nrejecttarget = nrejectout = 0;
    particle_type = box->forcefield->label2type[particle_label];
    label = "AVBMCOut";
}


/* ----------------------------------------------------------------------------
   Remove a random particle from the bonded region of another particle.
---------------------------------------------------------------------------- */

void AVBMCOut::perform_move()
{
    unsigned int i;

    //du = 0.;
    move_performed = false;

    // cannot remove particle if it does not exist
    if (box->npartype[particle_type] < 1) {
        nrejecttarget ++;
        return;
    }

    // pick target particle and count number of neighbors
    i = box->typeidx[particle_type][rng->next_int(box->npartype[particle_type])];
    std::vector<unsigned int> neigh_listi = box->distance_manager->build_neigh_list(i, r_abovesq);
    n_in = neigh_listi.size();

    // target atom needs neighbors in order to remove neighbor
    if (n_in < 1) {
        nrejectout ++;
        return;
    }

    // pick particle to be removed randomly among neighbors
    for (unsigned int j : rng->shuffle(neigh_listi)) {
        if (box->particles[j].type == particle_type) {
            du = -box->forcefield->comp_energy_par_force0(j); // this should be
            box->poteng += du;                                // baked into rm_particle
            particle_out = box->particles[j];                 // store old particle in
            box->rm_particle(j);                              // case move is rejected
            move_performed = true;
            return; //break;
        }
    }
}


/* ----------------------------------------------------------------------------
   Return the acceptance probability of move, given temperature 'temp' and
   chemical potential 'chempot'.
---------------------------------------------------------------------------- */

double AVBMCOut::accept(double temp, double chempot)
{
    if (!move_performed) return 0.;

    double dw = system->sampler->w(box->npar) - system->sampler->w(box->npar+1);
    double prefactor = n_in * box->npar / (v_in * (box->npar - 1));
    return prefactor * std::exp(-(du+chempot+dw)/temp);
}


/* ----------------------------------------------------------------------------
   Set back to old state if move is rejected
---------------------------------------------------------------------------- */

void AVBMCOut::reset()
{
    box->add_particle(particle_out);
    box->poteng -= du;
}


/* ----------------------------------------------------------------------------
   Update number of time this system size has occured if move was accepted
---------------------------------------------------------------------------- */

void AVBMCOut::update_size_histogram()
{
    box->update_size_histogram();
}


/* ----------------------------------------------------------------------------
   Represent move in a clean way
---------------------------------------------------------------------------- */

std::string AVBMCOut::repr()
{
    std::string move_info;
    move_info += "AVBMC particle deletion move\n";
    move_info += "    Radius of outer sphere: " + std::to_string(r_above) + "\n";
    move_info += "    Label of inserted atom: " + particle_label + "\n";
    return move_info;
}
