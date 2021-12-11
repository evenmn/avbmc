#include <iostream>
#include <cmath>
#include <vector>

#include "avbmcoutmol.h"
#include "../box.h"
#include "../particle.h"
#include "../molecule.h"



/* -----------------------------------------------------
   Aggregate-volume-biased out move. This is a 
   fictitous inter-box move, and works in the grand
   canonical ensemble only
-------------------------------------------------------- */

AVBMCOutMol::AVBMCOutMol(Box* box_in, const double r_above_in)
    : Moves(box_in)
{
    r_above = r_above_in;
    r_above2 = r_above * r_above;
    v_in = 1.; // 4 * pi * std::pow(r_above, 3)/3; // can be set to 1 according to Henrik
    not_accept = false;
}


/* ------------------------------------------------------
   Remove a random molecule from the bonded region of
   another similar molecule.
--------------------------------------------------------- */

void AVBMCOutMol::perform_move()
{
    // pick molecule to be removed
    int mol_idx = box->rng->choice(box->molecule_types->molecule_probs);
    molecule_out = box->molecule_types->construct_molecule(mol_idx);
    if(box->npar < 2 * molecule_out->natom){
        not_accept = true;
    }
    else{
        not_accept = false;
        du = -box->forcefield->comp_energy_mol(box->particles, molecule_out);
        particles_old = box->particles;
        for(int atom : molecule_out->atoms_idx){
            box->particles.erase(box->particles.begin() + atom);
        }
        box->npar -= molecule_out->natom;
    }
}


/* -------------------------------------------------------------
   Return the acceptance probability of move, given temperature
   'temp' and chemical potential 'chempot'.
---------------------------------------------------------------- */

double AVBMCOutMol::accept(double temp, double chempot)
{
    if(not_accept){
        return 0;
    }
    else{
        double dw = box->sampler->w(box->npar) - box->sampler->w(box->npar + molecule_out->natom);
        return n_in * box->npar / (v_in * (box->npar - molecule_out->natom)) * std::exp(-(du+chempot+dw)/temp);
    }
}


/* ------------------------------------------------------------
   Set back to old state if move is rejected
--------------------------------------------------------------- */

void AVBMCOutMol::reset()
{
    box->npar += molecule_out->natom;
    box->poteng -= du;
    box->particles = particles_old;
}
