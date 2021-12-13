#include <iostream>
#include <cmath>
#include <vector>

#include "avbmcout.h"
#include "../box.h"
#include "../particle.h"
#include "../molecule.h"



/* -----------------------------------------------------
   Aggregate-volume-biased out move. This is a 
   fictitous inter-box move, and works in the grand
   canonical ensemble only
-------------------------------------------------------- */

AVBMCOut::AVBMCOut(Box* box_in, const double r_above_in)
    : Moves(box_in)
{
    r_above = r_above_in;
    r_above2 = r_above * r_above;
    v_in = 1.; // 4 * pi * std::pow(r_above, 3)/3; // can be set to 1 according to Henrik
}


/* ------------------------------------------------------
   Remove a random particle from the bonded region of
   another particle.
--------------------------------------------------------- */

void AVBMCOut::perform_move()
{
    if(box->npar < 2){
        reject_move = true;
    }
    else{
        reject_move = false;
        // create local neighbor list of particle i

        int i = box->rng->next_int(box->npar);
        std::vector<int> neigh_listi = box->forcefield->build_neigh_list(i, r_above2);
        n_in = neigh_listi.size();

        if(n_in > 0){
            // pick particle to be removed
            int neigh_idx = box->rng->next_int(n_in);
            int j = neigh_listi[neigh_idx];  // particle to be removed
            du = -box->forcefield->comp_energy_par(box->particles, j);
            box->poteng += du;
            particle_out = box->particles[j];
            box->particles.erase(box->particles.begin() + j);
            box->npar --;
        }
        else{
            reject_move = true;
        }
    }
}


/* -------------------------------------------------------------
   Return the acceptance probability of move, given temperature
   'temp' and chemical potential 'chempot'.
---------------------------------------------------------------- */

double AVBMCOut::accept(double temp, double chempot)
{
    if(reject_move){
        return 0.;
    }
    else{
        double dw = box->sampler->w(box->npar) - box->sampler->w(box->npar+1);
        return n_in * box->npar / (v_in * (box->npar - 1)) * std::exp(-(du+chempot+dw)/temp);
    }
}


/* ------------------------------------------------------------
   Set back to old state if move is rejected
--------------------------------------------------------------- */

void AVBMCOut::reset()
{
    if (!reject_move){
        box->npar += 1;
        box->poteng -= du;
        box->particles.push_back(particle_out);
    }
}
