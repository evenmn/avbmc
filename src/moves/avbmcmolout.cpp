/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.

  Author(s): Even M. Nordhagen
  Email(s): evenmn@mn.uio.no
  Date: 2022-06-03 (last changed 2022-08-22)
---------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------
  AVBMC molecule deletion moves. The AVBMC moves were first proposed by Chen
  and Siepmann, "A Novel Monte Carlo Algorithm for Simulating Strongly
  Associating Fluids:  Applications to Water, Hydrogen Fluoride, and Acetic
  Acid" (2000). Non-bonded molecule moves have not yet been reported in the
  literature
---------------------------------------------------------------------------- */

#include <iostream>
#include <cmath>
#include <vector>

#include "avbmcmolout.h"
#include "../box.h"
#include "../system.h"
#include "../particle.h"
#include "../rng/rng.h"
#include "../distance_manager.h"
#include "../sampler/sampler.h"
#include "../boundary/boundary.h"
#include "../constraint/constraint.h"
#include "../forcefield/forcefield.h"


/* ----------------------------------------------------------------------------
   Aggregate-volume-biased molecule out move. This is a fictitous inter-box
   move, and works in the grand canonical ensemble only. 'molecule_in'
   specifies what molecule to remove. 'r_above_in' is the maximum distance 
   between target molecule (main atom) and molecule to remove (main atom). 
   'r_max_inner_in' is the maximum distance between main atom and other atoms
   to be considered a molecule. 'energy_bias_in' is a boolean controlling
   whether or not energy biasing should be used. 'target_mol_in' is a boolean
   defining the target molecule as the entire molecule or main atom.
---------------------------------------------------------------------------- */

AVBMCMolOut::AVBMCMolOut(System* system_in, Box* box_in,
         std::vector<Particle> molecule_in, const double r_above_in,
         const double r_inner_in, const bool energy_bias_in,
         const bool target_mol_in)
    : Moves(system_in)
{
    box = box_in;
    r_above = r_above_in;
    r_abovesq = r_above * r_above;
    r_inner = r_inner_in;
    energy_bias = energy_bias_in;
    target_mol = target_mol_in;
    molecule = molecule_in;
    natom = molecule.size();
    natom_inv = 1. / natom;
    v_in = 1.; // 4 * pi * std::pow(r_above, 3)/3;
    cum_time = 0.;
    naccept = ndrawn = nrejecttarget = nrejectout = 0;
    label = "AVBMCMolOut";

    std::string center = molecule[0].label;

    neigh_id_above = box->distance_manager->add_cutoff(r_above, center, center);
    neigh_id_inner = box->distance_manager->add_cutoff(r_inner);

    for (Particle &particle : molecule) {
        particle.type = box->forcefield->label2type.at(particle.label);
    }
    center_type = molecule[0].type;
}


/* ----------------------------------------------------------------------------
   Remove a random molecule from the bonded region of another similar molecule.
---------------------------------------------------------------------------- */

void AVBMCMolOut::perform_move()
{
    du = 0.;
    detected_out = false;

    // cannot remove particle if it does not exist
    if (box->npartype[center_type] < 1) {
        nrejecttarget ++;
        return;
    }

    // loop through potential target particles
    for (unsigned int t : rng->shuffle(box->typeidx[center_type])) {
        // loop through potential deletion particles
        std::vector<unsigned int> target_neigh = box->distance_manager->neigh_lists[neigh_id_above][t];
        // --- begin detect molecule
        for (unsigned int d : rng->shuffle(target_neigh)) {
            std::vector<unsigned int> molecule_idx_out, delete_neigh;
            if (box->particles[d].type == center_type) {
                molecule_idx_out.push_back(d);
                delete_neigh = box->distance_manager->neigh_lists[neigh_id_inner][d];
                n_in = delete_neigh.size();
                for (unsigned int i=1; i<natom; i++) {
                    detected_out = false;
                    for (unsigned int j=0; j<delete_neigh.size(); j++) {
                        if (molecule[i].type == box->particles[delete_neigh[j]].type) {
                            molecule_idx_out.push_back(delete_neigh[j]);
                            delete_neigh.erase(delete_neigh.begin() + j);
                            detected_out = true;
                            break;
                        }
                    }
                    if (!detected_out) break;
                }
            }
            // --- end detect molecule

            if (!detected_out) continue;

            //check_neigh_recu(d, 0, molecule_idx_out, delete_neigh);

            // compute change of energy when removing molecule
            if (box->store_energy) {
                // std::accumulate
                for (unsigned int k : molecule_idx_out) {
                    du -= box->forcefield->poteng_vec[k];
                }
            }
            else {
                for (unsigned int k : molecule_idx_out) {
                    du -= box->forcefield->comp_energy_par_force0(k);
                }
            }

            // remove molecule
            molecule_out.clear();
            std::sort(molecule_idx_out.begin(), molecule_idx_out.end(), std::greater<unsigned int>()); // sort in descending order
            for (unsigned int k : molecule_idx_out) {
                molecule_out.push_back(box->particles[k]);
                box->rm_particle(k);
            }
            box->poteng += du;
            nmolavg = n_in * natom_inv;
            return;
        }
    }
}


/* ----------------------------------------------------------------------------
   Return the acceptance probability of move, given temperature
   'temp' and chemical potential 'chempot'.
---------------------------------------------------------------------------- */

double AVBMCMolOut::accept(double temp, double chempot)
{
    if (!detected_out) {
        nrejectout ++;
        return 0.;
    }

    for (Constraint* constraint : box->constraints) {
        if (!constraint->verify()) return 0.;
    }

    double dw = system->sampler->w(box->npar) - system->sampler->w(box->npar + natom);
    double prefactor = nmolavg * box->npar / (v_in * (box->npar - natom)); 
    return prefactor * std::exp(-(du+chempot+dw)/temp);
}


/* ----------------------------------------------------------------------------
   Set back to old state if move is rejected
---------------------------------------------------------------------------- */

void AVBMCMolOut::reset()
{
    if (detected_out) {
        for (Particle particle : molecule_out) {
            box->add_particle(particle);
        }
        box->poteng -= du;
    }
}


/* ----------------------------------------------------------------------------
   Update number of time this system size has occured if move was accepted
---------------------------------------------------------------------------- */

void AVBMCMolOut::update_size_histogram()
{
    box->update_size_histogram();
}


/* ----------------------------------------------------------------------------
   Represent move in a clean way
---------------------------------------------------------------------------- */

std::string AVBMCMolOut::repr()
{
    std::string move_info;
    move_info += "AVBMC molecule deletion move\n";
    move_info += "    Radius of outer sphere: " + std::to_string(r_above) + "\n";
    move_info += "    Energy bias:            " + std::to_string(energy_bias) + "\n";
    move_info += "    Search target molecule: " + std::to_string(target_mol) + "\n";
    return move_info;
}
