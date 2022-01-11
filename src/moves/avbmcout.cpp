#include <iostream>
#include <cmath>
#include <vector>

#include "avbmcout.h"
#include "../box.h"
#include "../system.h"
#include "../particle.h"
#include "../molecule.h"
#include "../rng/rng.h"
#include "../sampler/sampler.h"
#include "../boundary/boundary.h"
#include "../forcefield/forcefield.h"


/* -----------------------------------------------------
   Aggregate-volume-biased out move. This is a 
   fictitous inter-box move, and works in the grand
   canonical ensemble only
-------------------------------------------------------- */

AVBMCOut::AVBMCOut(System* system_in, Box* box_in, std::string label_in,
                   const double r_above_in)
    : Moves(system_in)
{
    box = box_in;
    r_above = r_above_in;
    r_abovesq = r_above * r_above;
    v_in = 1.; // 4 * pi * std::pow(r_above, 3)/3;
    particle_label = label_in;
    label = "AVBMCOut";
}


/* ------------------------------------------------------
   Remove a random particle from the bonded region of
   another particle.
--------------------------------------------------------- */

void AVBMCOut::perform_move()
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


/* -------------------------------------------------------------
   Return the acceptance probability of move, given temperature
   'temp' and chemical potential 'chempot'.
---------------------------------------------------------------- */

double AVBMCOut::accept(double temp, double chempot)
{
    if (reject_move) {
        return 0.;
    }
    else {
        bool accept_boundary = box->boundary->correct_position();
        double dw = system->sampler->w(box->npar) - system->sampler->w(box->npar+1);
        return n_in * box->npar / (v_in * (box->npar - 1)) * std::exp(-(du+chempot+dw)/temp) * accept_boundary;
    }
}


/* ------------------------------------------------------------
   Set back to old state if move is rejected
--------------------------------------------------------------- */

void AVBMCOut::reset()
{
    if (!reject_move) {
        box->npar ++;
        box->poteng -= du;
        box->particles.push_back(particle_out);
    }
}


/* ----------------------------------------------------------
   Update number of time this system size has occured if
   move was accepted
------------------------------------------------------------- */

void AVBMCOut::update_nsystemsize()
{
    if (box->npar - 1 > box->nsystemsize.size()) {
        box->nsystemsize.resize(box->npar + 1);
    }
    box->nsystemsize[box->npar] ++;
}


/* -----------------------------------------------------
   Represent move in a clean way
-------------------------------------------------------- */

std::string AVBMCOut::repr()
{
    std::string move_info;
    move_info += "AVBMC particle deletion move\n";
    move_info += "    Radius of outer sphere: " + std::to_string(r_above) + "\n";
    return move_info;
}
