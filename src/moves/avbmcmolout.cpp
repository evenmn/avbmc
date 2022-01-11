#include <iostream>
#include <cmath>
#include <vector>

#include "avbmcmolout.h"
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

AVBMCMolOut::AVBMCMolOut(System* system_in, Box* box_in, const double r_above_in)
    : Moves(system_in)
{
    box = box_in;
    r_above = r_above_in;
    r_abovesq = r_above * r_above;
    v_in = 1.; // 4 * pi * std::pow(r_above, 3)/3; // can be set to 1 according to Henrik
    label = "AVBMCMolOut";
}


AVBMCMolOut::AVBMCMolOut(System* system_in, Box* box_in, std::vector<std::string> molecule_elements, const double r_above_in, const double r_max_inner_in)
    : Moves(system_in)
{
    box = box_in;
    r_above = r_above_in;
    r_abovesq = r_above * r_above;
    v_in = 1.; // 4 * pi * std::pow(r_above, 3)/3; // can be set to 1 according to Henrik
    label = "AVBMCMolOut";
}


/* ------------------------------------------------------
   Remove a random molecule from the bonded region of
   another similar molecule.
--------------------------------------------------------- */

void AVBMCMolOut::perform_move()
{
    reject_move = true;
    // pick molecule type to be removed
    MoleculeTypes* molecule_types = system->molecule_types;
    int mol_idx = rng->choice(molecule_types->molecule_probs);
    int natom = molecule_types->molecule_types[mol_idx].size();
    if (box->npar > 2 * natom){
        bool detected_out = false;
        int count = 0;
        while (!detected_out && count < box->npar){
            bool detected_target = false;
            std::vector<int> target_molecule = molecule_types->detect_molecule(box->particles, mol_idx, detected_target);
            if (detected_target){
                int i = target_molecule[0];  // center atom
                std::vector<int> neigh_listi = box->build_neigh_list(i, r_abovesq);  // neigh list of i
                int n_in = neigh_listi.size();  // number of neighbors of i
                if (n_in >= natom) {  // check that i has enough neighbors to potentially remove molecule
                    std::vector<Particle> particles;
                    for (int j=0; j < n_in; j++){
                        particles.push_back(box->particles[neigh_listi[j]]);
                    }
                    std::vector<int> molecule_out = molecule_types->detect_molecule(particles, mol_idx, detected_out);
                    if (detected_out) {
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
                        box->npar -= molecule_out.size();
                        nmolavg = (double) n_in / natom;
                        reject_move = false;
                        break;
                    }  // end 
                }  // end if n_in
            }  // end if constructed_target
            count ++;
        }  // end while
    }  // end if npar
    delete molecule_types;
}  // end perform_move


/* -------------------------------------------------------------
   Return the acceptance probability of move, given temperature
   'temp' and chemical potential 'chempot'.
---------------------------------------------------------------- */

double AVBMCMolOut::accept(double temp, double chempot)
{
    if (reject_move){
        return 0.;
    }
    else {
        double dw = system->sampler->w(box->npar) - system->sampler->w(box->npar + molecule_out->natom);
        return nmolavg * box->npar / (v_in * (box->npar - molecule_out->natom)) * std::exp(-(du+chempot+dw)/temp);
    }
}


/* ------------------------------------------------------------
   Set back to old state if move is rejected
--------------------------------------------------------------- */

void AVBMCMolOut::reset()
{
    if (!reject_move) {
        box->npar += molecule_out->natom;
        box->poteng -= du;
        box->particles = particles_old;
    }
}


/* ----------------------------------------------------------
   Update number of time this system size has occured if
   move was accepted
------------------------------------------------------------- */

void AVBMCMolOut::update_nsystemsize()
{
    if (box->npar + 1 > box->nsystemsize.size()) {
        box->nsystemsize.resize(box->npar + 1);
    }
    box->nsystemsize[box->npar] ++;
}


/* -----------------------------------------------------
   Represent move in a clean way
-------------------------------------------------------- */

std::string AVBMCMolOut::repr()
{
    std::string move_info;
    move_info += "AVBMC molecule deletion move\n";
    move_info += "    Radius of outer sphere: " + std::to_string(r_above) + "\n";
    return move_info;
}
