#include <iostream>
#include <cmath>
#include <vector>
#include <memory>

#include "avbmcout.h"
#include "../box.h"
#include "../system.h"
#include "../particle.h"
#include "../molecule.h"
#include "../rng/rng.h"
#include "../sampler/sampler.h"
#include "../forcefield/forcefield.h"


/* -----------------------------------------------------
   Aggregate-volume-biased out move. This is a 
   fictitous inter-box move, and works in the grand
   canonical ensemble only
-------------------------------------------------------- */

AVBMCOut::AVBMCOut(System* system_in, Box* box_in, const double r_above_in)
    : Moves(system_in)
{
    box = box_in;
    //boxes.push_back(box_in);
    r_above = r_above_in;
    r_abovesq = r_above * r_above;
    v_in = 1.; // 4 * pi * std::pow(r_above, 3)/3; // can be set to 1 according to Henrik
    label = "AVBMCOut";
}

/*
AVBMCOut::AVBMCOut(System* system_in, std::shared_ptr<Box> box_in, const double r_above_in)
    : Moves(system_in)
{
    box = box_in;
    boxes.push_back(box_in);
    r_above = r_above_in;
    r_abovesq = r_above * r_above;
    v_in = 1.; // 4 * pi * std::pow(r_above, 3)/3; // can be set to 1 according to Henrik
    label = "AVBMCOut";
}
*/

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
        int i = rng->next_int(box->npar);
        std::vector<int> neigh_listi = box->build_neigh_list(i, r_abovesq);
        n_in = neigh_listi.size();

        if(n_in > 0){
            // pick particle to be removed
            int neigh_idx = rng->next_int(n_in);
            int j = neigh_listi[neigh_idx];  // particle to be removed
            du = -system->forcefield->comp_energy_par(box->particles, j);
            box->poteng += du;
            particle_out = &box->particles[j];
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
        double dw = system->sampler->w(box->npar) - system->sampler->w(box->npar+1);
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
        box->particles.emplace_back(*particle_out);
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
