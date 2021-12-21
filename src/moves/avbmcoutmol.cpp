#include <iostream>
#include <cmath>
#include <vector>

#include "avbmcoutmol.h"
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

AVBMCOutMol::AVBMCOutMol(System* system_in, Box* box_in, const double r_above_in)
    : Moves(system_in)
{
    box = box_in;
    boxes.push_back(box_in);
    r_above = r_above_in;
    r_abovesq = r_above * r_above;
    v_in = 1.; // 4 * pi * std::pow(r_above, 3)/3; // can be set to 1 according to Henrik
    label = "AVBMCOutMol";
}


/* ------------------------------------------------------
   Remove a random molecule from the bonded region of
   another similar molecule.
--------------------------------------------------------- */

void AVBMCOutMol::perform_move()
{
    reject_move = true;
    // pick molecule type to be removed
    MoleculeTypes* molecule_types = system->molecule_types;
    int mol_idx = rng->choice(molecule_types->molecule_probs);
    int natom = molecule_types->molecule_types[mol_idx].size();
    if (box->npar > 2 * natom){
        bool constructed_out = false;
        int count = 0;
        while (!constructed_out && count < box->npar){
            bool constructed_target = false;
            Molecule* molecule_target = molecule_types->construct_molecule(box->particles, mol_idx, constructed_target);
            if (constructed_target){
                int i = molecule_target->atoms_idx[0];  // center atom
                std::vector<int> neigh_listi = box->build_neigh_list(i, r_abovesq);  // neigh list of i
                int n_in = neigh_listi.size();  // number of neighbors of i

                if (n_in >= natom) {  // check that i has enough neighbors to potentially remove molecule
                    std::vector<Particle *> particles;
                    for (int j=0; j < n_in; j++){
                        particles.push_back(box->particles[neigh_listi[j]]);
                    }
                    molecule_out = molecule_types->construct_molecule(particles, mol_idx, constructed_out);
                    if (constructed_out) {
                        du = -system->forcefield->comp_energy_mol(box->particles, molecule_out);
                        particles_old = box->particles;
                        std::vector<int> atoms_idx = molecule_out->atoms_idx;
                        std::sort(atoms_idx.begin(), atoms_idx.end(), std::greater<int>()); // sort in decending order
                        for(int atom_idx : atoms_idx){
                            box->particles.erase(box->particles.begin() + atom_idx);
                        }
                        box->poteng += du;
                        box->npar -= molecule_out->natom;
                        nmolavg = (double) n_in / natom;
                        reject_move = false;
                        break;
                    }  // end 
                }  // end if n_in
            }  // end if constructed_target
            count ++;
        }  // end while
    }  // end if npar
}  // end perform_move


/* -------------------------------------------------------------
   Return the acceptance probability of move, given temperature
   'temp' and chemical potential 'chempot'.
---------------------------------------------------------------- */

double AVBMCOutMol::accept(double temp, double chempot)
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

void AVBMCOutMol::reset()
{
    if (!reject_move) {
        box->npar += molecule_out->natom;
        box->poteng -= du;
        box->particles = particles_old;
    }
}


/* -----------------------------------------------------
   Represent move in a clean way
-------------------------------------------------------- */

std::string AVBMCOutMol::repr()
{
    std::string move_info;
    move_info += "AVBMC molecule deletion move\n";
    move_info += "    Radius of outer sphere: " + std::to_string(r_above) + "\n";
    return move_info;
}
