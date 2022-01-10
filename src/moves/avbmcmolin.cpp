/* Information */

/* -----------------------------------------------------
   Aggregate-volume-biased in move. This is a fictitious 
   inter-box move, and works in the grand canonical 
   essemble only. A predefined molecule is inserted
-------------------------------------------------------- */

#include <iostream>
#include <cmath>
#include <valarray>
#include <vector>

#include "avbmcmolin.h"
#include "../box.h"
#include "../system.h"
#include "../particle.h"
#include "../molecule.h"
#include "../rng/rng.h"
#include "../sampler/sampler.h"
#include "../forcefield/forcefield.h"


/* -----------------------------------------------------
   AVBMC molecule insertation move constructor. Molecule
   is inserted in the bonded region of an equivalent
   molecule, defined by inner radius 'r_below_in' and
   outer radius 'r_above_in'. 
-------------------------------------------------------- */

AVBMCMolIn::AVBMCMolIn(System* system_in, Box* box_in,
                       const double r_below_in, const double r_above_in)
    : Moves(system_in)
{
    box = box_in;
    r_below = r_below_in;
    r_above = r_above_in;
    r_abovesq = r_above * r_above;
    r_belowsq = r_below * r_below;
    v_in = 1.; // 4 * pi * std::pow(r_above, 3)/3;
    label = "AVBMCMolIn";
}


/* ------------------------------------------------------------
   Insert molecule into the bonded region of an equivalent
   molecule
--------------------------------------------------------------- */

void AVBMCMolIn::perform_move()
{
    // pick molecule type to be inserted
    MoleculeTypes* molecule_types = system->molecule_types;
    std::cout << "LENGTH: " << system->molecule_types->molecule_probs.size() << std::endl;
    for (double prob : molecule_types->molecule_probs) {
        std::cout << "prob: " << prob << std::endl;
    }
    int mol_idx = rng->choice(molecule_types->molecule_probs);
    std::vector<std::valarray<double> > positions = molecule_types->default_mols[mol_idx];
    natom = molecule_types->molecule_elements[mol_idx].size();

    // rotate molecule arbitrary
    positions = rotate_molecule(positions);

    // detect target molecule
    bool detected = false;
    std::vector<int> target_molecule = molecule_types->detect_molecule(box->particles, mol_idx, detected);
    if (!detected) {
        reject_move = true;
        std::cout << "rejected" << std::endl;
    }
    else {
        reject_move = false;
        int i = target_molecule[0];  // center atom
        std::vector<int> neigh_listi = box->build_neigh_list(i, r_abovesq);
        int n_in = neigh_listi.size();
        nmolavg = (double) n_in / natom;

        // shift in-molecule relative to target molecule
        std::valarray<double> dr(system->ndim);
        double normsq = norm(dr);
        while (normsq > r_abovesq || normsq < r_belowsq) {
            for (double &d : dr) {
                d = r_above * (2 * rng->next_double() - 1);
            }
            normsq = norm(dr);
        }

        // construct new particles
        for (int j=0; j < natom; j++) {
            positions[j] += box->particles[i].r + dr;
            std::string element = molecule_types->molecule_elements[mol_idx][j];
            Particle particle(element, positions[j]);
            particle.type = system->label2type[element];
            box->particles.push_back(particle);
            box->npar ++;
        }

        // compute energy difference
        du = 0.;
        for (int j=0; j < natom; j++) {
            du += system->forcefield->comp_energy_par(box->particles, box->npar - j - 1);
        }
        box->poteng += du;
    }
    delete molecule_types;
    molecule_types = nullptr;
}


/* ---------------------------------------------------------
   Get acceptance probability of move
------------------------------------------------------------ */

double AVBMCMolIn::accept(double temp, double chempot)
{
    if (reject_move) {
        return 0.;
    }
    else {
        double dw = system->sampler->w(box->npar) - system->sampler->w(box->npar - natom);
        return (v_in * box->npar) / ((nmolavg + 1) * (box->npar + natom)) * std::exp(-(du-chempot+dw)/temp);
    }
}


/* ----------------------------------------------------------
   Set back to old state before move
------------------------------------------------------------- */

void AVBMCMolIn::reset()
{
    if (!reject_move) {
        box->npar -= natom;
        box->poteng -= du;
        box->particles.erase(box->particles.begin() + box->npar);
    }
}


/* ----------------------------------------------------------
   Update number of time this system size has occured if
   move was accepted
------------------------------------------------------------- */

void AVBMCMolIn::update_nsystemsize()
{
    if (box->npar + 1 > box->nsystemsize.size()) {
        box->nsystemsize.resize(box->npar + 1);
    }
    box->nsystemsize[box->npar] ++;
}


/* -----------------------------------------------------
   Represent move in a clean way
-------------------------------------------------------- */

std::string AVBMCMolIn::repr()
{
    std::string move_info;
    move_info += "AVBMC molecule insertion move\n";
    move_info += "    Radius of outer sphere: " + std::to_string(r_above) + "\n";
    move_info += "    Radius of inner sphere: " + std::to_string(r_below) + "\n";
    return move_info;
}
