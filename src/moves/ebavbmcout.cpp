/* ----------------------------------------------------------------------------
------------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------
   Energy-biased Aggregation-volume-biased Monte Carlo (EB-AVBMC) deletion
   move. This is a fictitious inter-box move, and works in the grand-canonical
   essemble only. To maintain detailed balance, this move always have to be used
   in conjunction with an EB-AVBMC insertion move, with the same radius of the 
   bonded region and the same move probability. This is ensured when applying
   the 'EB-AVBMC' move. The EB-AVBMC moves were first proposed by Loeffler, 
   Sepehri and Chen (2015).
------------------------------------------------------------------------------- */

#include <iostream>
#include <cmath>
#include <vector>

#include "ebavbmcout.h"
#include "../box.h"
#include "../system.h"
#include "../particle.h"
#include "../rng/rng.h"
#include "../sampler/sampler.h"
#include "../boundary/boundary.h"
#include "../forcefield/forcefield.h"


/* ----------------------------------------------------------------------------
   EB-AVBMCOut constructor
------------------------------------------------------------------------------- */

EBAVBMCOut::EBAVBMCOut(System* system_in, Box* box_in, std::string label_in,
                         const double r_above_in)
    : Moves(system_in)
{
    box = box_in;
    r_above = r_above_in;
    r_abovesq = r_above * r_above;
    v_in = 1.; // 4 * pi * std::pow(r_above, 3)/3;
    particle_label = label_in;
    label = "EB-AVBMCOut";
}


/* ----------------------------------------------------------------------------
   Remove a random particle from the bonded region of
   another particle.
------------------------------------------------------------------------------- */

void EBAVBMCOut::perform_move()
{
    unsigned int i, j, counti, countj, neigh_idx;
    reject_move = true;
    if (box->npar > 2) {
        counti = 0;
        // ensure that removed particle is of correct type
        while (reject_move && counti < box->npar) {
            i = rng->next_int(box->npar);
            if (box->particles[i].label == particle_label) {
                reject_move = false;
            }
            counti ++;
        }
        if (!reject_move) {
            std::vector<int> neigh_listi = box->build_neigh_list(i, r_abovesq);
            n_in = neigh_listi.size();

            if (n_in > 0) {
                // pick particle to be removed randomly among neighbors
                reject_move = true;
                countj = 0;
                // ensure that removed particle is of correct type
                while (reject_move && countj < n_in) {
                    neigh_idx = rng->next_int(n_in);
                    j = neigh_listi[neigh_idx];  // particle to be removed
                    if (box->particles[j].label == particle_label) {
                        reject_move = false;
                    }
                    countj ++;
                }
                if (!reject_move) {
                    du = -system->forcefield->comp_energy_par(box->particles, j);
                    box->poteng += du;
                    particle_out = box->particles[j];
                    box->particles.erase(box->particles.begin() + j);
                    box->npar --;
                }
            }
        }
    }
}


/* ----------------------------------------------------------------------------
   Return the acceptance probability of move, given temperature
   'temp' and chemical potential 'chempot'.
------------------------------------------------------------------------------- */

double EBAVBMCOut::accept(double temp, double chempot)
{
    if (reject_move) {
        return 0.;
    }
    else {
        bool accept_boundary = box->boundary->correct_position();
        double dw = system->sampler->w(box->npar) - system->sampler->w(box->npar+1);
        double prefactor = n_in * box->npar / (v_in * (box->npar - 1));
        return prefactor * std::exp(-(du+chempot+dw)/temp) * accept_boundary;
    }
}


/* ----------------------------------------------------------------------------
   Set back to old state if move is rejected
------------------------------------------------------------------------------- */

void EBAVBMCOut::reset()
{
    if (!reject_move) {
        box->npar ++;
        box->poteng -= du;
        box->particles.push_back(particle_out);
    }
}


/* ----------------------------------------------------------------------------
   Update number of time this system size has occured if
   move was accepted
------------------------------------------------------------------------------- */

void EBAVBMCOut::update_nsystemsize()
{
    if (box->npar - 1 > box->nsystemsize.size()) {
        box->nsystemsize.resize(box->npar + 1);
    }
    box->nsystemsize[box->npar] ++;
}


/* ----------------------------------------------------------------------------
   Represent move in a clean way
------------------------------------------------------------------------------- */

std::string EBAVBMCOut::repr()
{
    std::string move_info;
    move_info += "EB-AVBMC particle deletion move\n";
    move_info += "    Radius of outer sphere: " + std::to_string(r_above) + "\n";
    move_info += "    Label of inserted atom: " + particle_label + "\n";
    return move_info;
}
