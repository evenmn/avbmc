/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.

  Author(s): Even M. Nordhagen
  Email(s): evenmn@mn.uio.no
  Date: 2022-06-03 (last changed 2022-06-03)
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
    naccept = ndrawn = 0;
    label = "AVBMCMolOut";

    // neigh_id_above = box->distance_manager->add_cutoff(r_above,
    // particles[0].label, particles[0].label, false); neigh_id_inner =
    // box->distance_manager->add_cutoff(r_inner, particles[0].label,
    // particles[1].label, false);
    //neigh_id_above = box->distance_manager->add_cutoff(r_above);
    //neigh_id_inner = box->distance_manager->add_cutoff(r_inner);

    for (Particle &particle : molecule) {
        particle.type = box->forcefield->label2type.at(particle.label);
    }
}


/* ----------------------------------------------------------------------------
   Detect target molecule (or atom if target_mol=false). Make maximum 'npar' 
   attempts of finding target molecule.
---------------------------------------------------------------------------- */
/*
unsigned int AVBMCMolOut::detect_target_molecule(bool &detected)
{
    unsigned int count, i;
    std::vector<unsigned int> target_molecule;
    std::vector<std::vector<unsigned int> > neigh_list_inner;

    count = i = 0;
    while (count < box->npar || detected) {
        count ++;
        if (target_mol) {
            neigh_list_inner = box->distance_manager->neigh_lists[neigh_id_inner];
            target_molecule = detect_molecule(neigh_list_inner, molecule, detected);
            i = target_molecule[0];
        }
        else {
            i = rng->next_int(box->npar);
            detected = (box->particles[i].type == molecule[0].type);
        }
    }
    return i;
}
*/

/* ----------------------------------------------------------------------------
   Detect deletion molecule.
---------------------------------------------------------------------------- */

std::vector<unsigned int> AVBMCMolOut::detect_deletion_molecule(unsigned int i,
    bool &detected)
{
    std::vector<unsigned int> molecule_out, molecule_out2, neigh_listi;
    molecule_out.clear();

    neigh_listi = box->distance_manager->build_neigh_list(i, r_abovesq);
    //neigh_listi = box->distance_manager->neigh_lists[neigh_id_above][i];
    n_in = neigh_listi.size();
    if (n_in < natom) {  // ensure that there is a least one molecule left
        detected = false;
        return molecule_out2;
    }
    std::vector<Particle> particles_tmp;
    for (unsigned int j=0; j < n_in; j++){
        particles_tmp.push_back(box->particles[neigh_listi[j]]);
    }

    //neigh_list_inner = box->distance_manager->neigh_lists[neigh_id_inner];
    //molecule_out2 = detect_molecule(neigh_list_inner, particles_tmp, molecule, detected);
    molecule_out2 = detect_molecule(particles_tmp, molecule, detected, r_inner, box);
    // std::transform
    for (unsigned int idx : molecule_out2) {
        molecule_out.push_back(neigh_listi[idx]);
        //molecule_out2.push_back(particles_tmp[idx]);
    }
    return molecule_out;
}


/* ----------------------------------------------------------------------------
   Remove a random molecule from the bonded region of another similar molecule.
---------------------------------------------------------------------------- */

void AVBMCMolOut::perform_move()
{
    unsigned int i, n_in;
    std::vector<unsigned int> molecule_idx_out;

    detected_out = false;

    // cannot remove particle if it does not exist
    if (box->npartype[molecule[0].type] < 1) {
        return;
    }

    //if (box->npar < 2 * natom - 1) {
    //    // do not remove molecule if there is less than two molecules available
    //    return;
    //}

    // detect target particle (molecule)
    /*
    if (target_mol) {
        i = detect_target_molecule(detected_target);
        if (!detected_target) {
            return;
        }
    }
    else {
        i = box->typeidx[particles[0].type][rng->next_int(box->npartype[particles[0].type])];
    }
    */
    i = box->typeidx[molecule[0].type][rng->next_int(box->npartype[molecule[0].type])];
    std::vector<unsigned int> neigh_listi = box->distance_manager->build_neigh_list(i, r_abovesq);
    n_in = neigh_listi.size(); 
    nmolavg = n_in * natom_inv;

    // target particle needs at least natom neighbors
    if (n_in < natom) {
        return;
    }
    molecule_idx_out = detect_deletion_molecule(i, detected_out);
    if (!detected_out) return;

    // compute change of energy when removing molecule
    du = 0.;
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
}


/* ----------------------------------------------------------------------------
   Return the acceptance probability of move, given temperature
   'temp' and chemical potential 'chempot'.
---------------------------------------------------------------------------- */

double AVBMCMolOut::accept(double temp, double chempot)
{
    if (!detected_out) return 0.;

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
