#include <iostream>
#include <cmath>
#include <valarray>
#include <vector>
#include <memory>

#include "avbmcin.h"
#include "../box.h"
#include "../system.h"
#include "../particle.h"
#include "../molecule.h"
#include "../rng/rng.h"
#include "../sampler/sampler.h"
#include "../forcefield/forcefield.h"


/* -----------------------------------------------------
   Aggregate-volume-biased in move. This is a fictitious 
   inter-box move, and works in the grand canonical 
   essemble only. A defined molecule is inserted
-------------------------------------------------------- */


AVBMCIn::AVBMCIn(System* system_in, Box* box_in, std::string particle_label,
                 const double r_below_in, const double r_above_in)
    : Moves(system_in)
{
    box = box_in;
    //boxes.push_back(box_in);
    r_below = r_below_in;
    r_above = r_above_in;
    r_abovesq = r_above * r_above;
    r_belowsq = r_below * r_below;
    v_in = 1.; // 4 * pi * std::pow(r_above, 3)/3; // can be set to 1 according to Henrik

    //type = 0;      // type and label of inserted particle
    label_in = particle_label;
    label = "AVBMCIn ";
}
/*
AVBMCIn::AVBMCIn(System* system_in, std::shared_ptr<Box> box_in, const double r_below_in, const double r_above_in)
    : Moves(system_in)
{
    box = box_in;
    boxes.push_back(box_in);
    r_below = r_below_in;
    r_above = r_above_in;
    r_abovesq = r_above * r_above;
    r_belowsq = r_below * r_below;
    v_in = 1.; // 4 * pi * std::pow(r_above, 3)/3; // can be set to 1 according to Henrik

    //type = 0;      // type and label of inserted particle
    label_in = "Ar";  // this has to be generalized
    label = "AVBMCIn ";
}


AVBMCIn::AVBMCIn(System* system_in, std::shared_ptr<Box> box_in, std::string particle_label,
                 const double r_below_in, const double r_above_in)
    : Moves(system_in)
{
    box = box_in;
    boxes.push_back(box_in);
    r_below = r_below_in;
    r_above = r_above_in;
    r_abovesq = r_above * r_above;
    r_belowsq = r_below * r_below;
    v_in = 1.; // 4 * pi * std::pow(r_above, 3)/3; // can be set to 1 according to Henrik

    label_in = particle_label;
    label = "AVBMCIn ";
}*/
/*
AVBMCIn::AVBMCIn(System* system_in, Box* box_in, const double r_below_in, const double r_above_in)
    : Moves(system_in)
{
    box = box_in;
    boxes.push_back(box_in);
    r_below = r_below_in;
    r_above = r_above_in;
    r_abovesq = r_above * r_above;
    r_belowsq = r_below * r_below;
    v_in = 1.; // 4 * pi * std::pow(r_above, 3)/3; // can be set to 1 according to Henrik

    type = 0;      // type and label of inserted particle
    label_in = "Ar";  // this has to be generalized
    label = "AVBMCIn ";
}
*/

/* ------------------------------------------------------------
   Insert molecule into the bonded region 
--------------------------------------------------------------- */

void AVBMCIn::perform_move()
{
    // pick particle i and create local neighbor list of particle
    int i = rng->next_int(box->npar);
    std::vector<int> neigh_listi = box->build_neigh_list(i, r_abovesq);
    n_in = neigh_listi.size();

    // construct new particle
    std::valarray<double> dr(system->ndim);
    double normsq = norm(dr);
    while(normsq > r_abovesq || normsq < r_belowsq){
        for(double &d : dr){
            d = 2 * rng->next_double() - 1;
        }
        normsq = norm(dr);
    }

    Particle particle_in = Particle(label_in, box->particles[i].r + dr);
    //box->npar ++;
    particle_in.type = system->label2type.at(label_in);
    box->particles.push_back(particle_in);

    // compute du
    du = system->forcefield->comp_energy_par(box->particles, box->npar - 1);
    box->poteng += du;
}


/* ---------------------------------------------------------
   Get acceptance probability of move
------------------------------------------------------------ */

double AVBMCIn::accept(double temp, double chempot)
{
    double dw = system->sampler->w(box->npar) - system->sampler->w(box->npar-1);
    return (v_in * box->npar) / ((n_in + 1) * (box->npar + 1)) * std::exp(-(du-chempot+dw)/temp);
}


/* ----------------------------------------------------------
   Set back to old state before move is move was rejected
------------------------------------------------------------- */

void AVBMCIn::reset()
{
    box->npar --;
    box->poteng -= du;
    box->particles.erase(box->particles.begin() + box->npar);
}


/* ----------------------------------------------------------
   Update number of time this system size has occured if
   move was accepted
------------------------------------------------------------- */

void AVBMCIn::update_nsystemsize()
{
    if (box->npar - 1 > box->nsystemsize.size()) {
        box->nsystemsize.resize(box->npar + 1);
    }
    box->nsystemsize[box->npar] ++;
}


/* -----------------------------------------------------
   Represent move in a clean way
-------------------------------------------------------- */

std::string AVBMCIn::repr()
{
    std::string move_info;
    move_info += "AVBMC particle insertion move\n";
    move_info += "    Radius of outer sphere: " + std::to_string(r_above) + "\n";
    move_info += "    Radius of inner sphere: " + std::to_string(r_below) + "\n";
    return move_info;
}
