#include <iostream>
#include <cmath>
#include <valarray>
#include <vector>

#include "avbmcin.h"
#include "../box.h"
#include "../particle.h"
#include "../molecule.h"


/* -----------------------------------------------------
   Aggregate-volume-biased in move. This is a fictitious 
   inter-box move, and works in the grand canonical 
   essemble only. A defined molecule is inserted
-------------------------------------------------------- */

AVBMCIn::AVBMCIn(Box* box_in, const double r_below_in, const double r_above_in)
    : Moves(box_in)
{
    r_below = r_below_in;
    r_above = r_above_in;
    r_abovesq = r_above * r_above;
    r_belowsq = r_below * r_below;
    v_in = 1.; // 4 * pi * std::pow(r_above, 3)/3; // can be set to 1 according to Henrik

    type = 0;      // type and label of inserted particle
    label = "Ar";  // this has to be generalized
}


/* ------------------------------------------------------------
   Insert molecule into the bonded region 
--------------------------------------------------------------- */

void AVBMCIn::perform_move()
{
    // pick particle i and create local neighbor list of particle
    int i = box->rng->next_int(box->npar);
    std::vector<int> neigh_listi = box->forcefield->build_neigh_list(i, r_abovesq);
    n_in = neigh_listi.size();

    /*
    // compute norm
    auto norm = [] (std::valarray<double> x) -> double { 
        double sqrd_sum = 0.;
        for(double x_ : x){
            sqrd_sum += x_ * x_;
        }
        return sqrd_sum;
    };
    */

    // construct new particle
    std::valarray<double> dr(box->ndim);
    double normsq = std::pow(dr, 2).sum(); //norm(dr);
    while(normsq > r_abovesq || normsq < r_belowsq){
        for(double &d : dr){
            d = 2 * box->rng->next_double() - 1;
        }
        normsq = std::pow(dr, 2).sum();
    }

    Particle *particle_in = new Particle(label, box->particles[i]->r + dr);
    particle_in->type = type;
    box->particles.push_back(particle_in);
    box->npar ++;

    // compute du
    du = box->forcefield->comp_energy_par(box->particles, box->npar - 1);
    box->poteng += du;
}


/* ---------------------------------------------------------
   Get acceptance probability of move
------------------------------------------------------------ */

double AVBMCIn::accept(double temp, double chempot)
{
    double dw = box->sampler->w(box->npar) - box->sampler->w(box->npar-1);
    return (v_in * box->npar) / ((n_in + 1) * (box->npar + 1)) * std::exp(-(du-chempot+dw)/temp);
}


/* ----------------------------------------------------------
   Set back to old state before move
------------------------------------------------------------- */

void AVBMCIn::reset()
{
    box->npar --;
    box->poteng -= du;
    box->particles.erase(box->particles.begin() + box->npar);
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
