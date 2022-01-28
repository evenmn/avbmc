#include <iostream>
#include <cmath>
#include <vector>

#include "avbmcmoloutres.h"
#include "../box.h"
#include "../system.h"
#include "../particle.h"
#include "../rng/rng.h"
#include "../sampler/sampler.h"
#include "../boundary/boundary.h"
#include "../forcefield/forcefield.h"


/* ----------------------------------------------------------------------------
   Aggregate-volume-biased out move. This is a fictitous inter-box move, and
   works in the grand canonical ensemble only
------------------------------------------------------------------------------- */

AVBMCMolOutRes::AVBMCMolOutRes(System* system_in, Box* box_in,
         std::vector<Particle> molecule_in, const double r_above_in,
         const double r_max_inner_in)
    : Moves(system_in)
{
    box = box_in;
    r_above = r_above_in;
    r_abovesq = r_above * r_above;
    r_max_inner = r_max_inner_in;
    molecule = molecule_in;
    natom = molecule.size();
    natom_inv = 1. / natom;
    v_in = 1.; // 4 * pi * std::pow(r_above, 3)/3; // can be set to 1 according to Henrik
    label = "AVBMCMolOutRes";
}


/* ----------------------------------------------------------------------------
   Remove a random molecule from the bonded region of another similar molecule.
------------------------------------------------------------------------------- */

void AVBMCMolOutRes::perform_move()
{
    bool detected_out;
    unsigned int count, i, n_in;
    reject_move = true;
    if (box->npar > 2 * natom - 1) {
        count = 0;
        while (count < box->npar && reject_move) {
            count ++;
            i = rng->next_int(box->npar);
            if (box->particles[i].type == molecule[0].type) {
                std::vector<int> neigh_listi = box->build_neigh_list(i, r_abovesq);  // neigh list of i
                n_in = neigh_listi.size();  // number of neighbors of i
                if (n_in >= natom) {  // check that i has enough neighbors to potentially remove molecule
                    std::vector<Particle> particles;
                    for (unsigned int j=0; j < n_in; j++){
                        particles.push_back(box->particles[neigh_listi[j]]);
                    }
                    detected_out = false;
                    std::vector<int> molecule_out = detect_molecule(particles, molecule, detected_out, r_max_inner);
                    if (detected_out) {
                        reject_move = false;
                        // compute change of energy when removing molecule
                        du = 0.;
                        for (int j : molecule_out) {
                            du -= system->forcefield->comp_energy_par(box->particles, i);
                        }
                        // remove molecule
                        particles_old = box->particles;
                        std::sort(molecule_out.begin(), molecule_out.end(), std::greater<int>()); // sort in decending order
                        for (int j : molecule_out){
                            box->particles.erase(box->particles.begin() + j);
                        }
                        box->poteng += du;
                        box->npar -= natom;
                        nmolavg = n_in * natom_inv;
                    }  // end 
                }  // end if n_in
            }  // end if constructed_target
        }  // end while
    }  // end if npar
}  // end perform_move


/* ----------------------------------------------------------------------------
   Return the acceptance probability of move, given temperature
   'temp' and chemical potential 'chempot'.
------------------------------------------------------------------------------- */

double AVBMCMolOutRes::accept(double temp, double chempot)
{
    if (reject_move){
        return 0.;
    }
    else {
        bool accept_boundary = box->boundary->correct_position();
        double dw = system->sampler->w(box->npar) - system->sampler->w(box->npar + natom);
        double prefactor = nmolavg * box->npar / (v_in * (box->npar - natom)); 
        return prefactor * accept_boundary * std::exp(-(du+chempot+dw)/temp);
    }
}


/* ----------------------------------------------------------------------------
   Set back to old state if move is rejected
------------------------------------------------------------------------------- */

void AVBMCMolOutRes::reset()
{
    if (!reject_move) {
        box->npar += natom;
        box->poteng -= du;
        box->particles = particles_old;
    }
}


/* ----------------------------------------------------------------------------
   Update number of time this system size has occured if
   move was accepted
------------------------------------------------------------------------------- */

void AVBMCMolOutRes::update_nsystemsize()
{
    if (box->npar + 1 > box->nsystemsize.size()) {
        box->nsystemsize.resize(box->npar + 1);
    }
    box->nsystemsize[box->npar] ++;
}


/* ----------------------------------------------------------------------------
   Represent move in a clean way
------------------------------------------------------------------------------- */

std::string AVBMCMolOutRes::repr()
{
    std::string move_info;
    move_info += "AVBMC molecule deletion move\n";
    move_info += "    Radius of outer sphere: " + std::to_string(r_above) + "\n";
    return move_info;
}
