/* ----------------------------------------------------------------------------
------------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------
   Energy-biased Aggregation-volume-biased Monte Carlo (EB-AVBMC) insertion
   move. This is a fictitious inter-box move, and works in the grand-canonical
   essemble only. To maintain detailed balance, this move always have to be
   used in conjunction with an EB-AVBMC deletion move, with the same radius of
   the bonded region and the same move probability. This is ensured when
   applying the 'EB-AVBMC' move. The EB-AVBMC moves were first proposed by
   Toeffler, Sepehri and Chen (2015).
------------------------------------------------------------------------------- */

#include <iostream>
#include <cmath>
#include <valarray>
#include <vector>

#include "ebavbmcin.h"
#include "../box.h"
#include "../system.h"
#include "../particle.h"
#include "../rng/rng.h"
#include "../sampler/sampler.h"
#include "../boundary/boundary.h"
#include "../forcefield/forcefield.h"


/* ----------------------------------------------------------------------------
   EB-AVBMCIn constructor
------------------------------------------------------------------------------- */


EBAVBMCIn::EBAVBMCIn(System* system_in, Box* box_in, std::string label_in,
                 const double r_below_in, const double r_above_in)
    : Moves(system_in)
{
    box = box_in;
    r_below = r_below_in;
    r_above = r_above_in;
    r_abovesq = r_above * r_above;
    r_belowsq = r_below * r_below;
    v_in = 1.; // 4 * pi * std::pow(r_above, 3)/3; // can be set to 1 according to Henrik

    particle_label = label_in;
    particle_type = system->forcefield->label2type.at(label_in);
    label = "EB-AVBMCIn ";
}


/* ----------------------------------------------------------------------------
   Insert molecule into the bonded region 
------------------------------------------------------------------------------- */

void EBAVBMCIn::perform_move()
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

    Particle particle_in(particle_label, box->particles[i].r + dr);
    box->npar ++;
    particle_in.type = particle_type;
    box->particles.push_back(particle_in);

    // compute du
    du = system->forcefield->comp_energy_par(box->particles, box->npar - 1);
    box->poteng += du;
}


/* ----------------------------------------------------------------------------
   Get acceptance probability of move
------------------------------------------------------------------------------- */

double EBAVBMCIn::accept(double temp, double chempot)
{
    bool accept_boundary = box->boundary->correct_position();
    double dw = system->sampler->w(box->npar) - system->sampler->w(box->npar-1);
    double prefactor = v_in * box->npar / ((n_in + 1) * (box->npar + 1));
    return prefactor * std::exp(-(du-chempot+dw)/temp) * accept_boundary;
}


/* ----------------------------------------------------------------------------
   Set back to old state before move is move was rejected
------------------------------------------------------------------------------- */

void EBAVBMCIn::reset()
{
    box->npar --;
    box->poteng -= du;
    box->particles.erase(box->particles.begin() + box->npar);
}


/* ----------------------------------------------------------------------------
   Update number of time this system size has occured if move was accepted
------------------------------------------------------------------------------- */

void EBAVBMCIn::update_nsystemsize()
{
    if (box->npar + 1 > box->nsystemsize.size()) {
        box->nsystemsize.resize(box->npar + 1);
    }
    box->nsystemsize[box->npar] ++;
}


/* ----------------------------------------------------------------------------
   Represent move in a clean way
------------------------------------------------------------------------------- */

std::string EBAVBMCIn::repr()
{
    std::string move_info;
    move_info += "EB-AVBMC particle insertion move\n";
    move_info += "    Radius of outer sphere: " + std::to_string(r_above) + "\n";
    move_info += "    Radius of inner sphere: " + std::to_string(r_below) + "\n";
    move_info += "    Label of inserted atom: " + particle_label + "\n";
    return move_info;
}