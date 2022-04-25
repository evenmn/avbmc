/* Information */

/* ----------------------------------------------------------------------------
   Aggregate-volume-biased in move. This is a fictitious inter-box move, and
   works in the grand canonical essemble only. A predefined molecule is
   inserted.
------------------------------------------------------------------------------- */

#include <iostream>
#include <cmath>
#include <valarray>
#include <vector>

#include "avbmcmolin.h"
#include "../box.h"
#include "../system.h"
#include "../particle.h"
#include "../rng/rng.h"
#include "../sampler/sampler.h"
#include "../boundary/boundary.h"
#include "../forcefield/forcefield.h"


/* ----------------------------------------------------------------------------
   AVBMC molecule insertation move constructor. Molecule is inserted in the
   bonded region of an equivalent molecule, defined by inner radius
   'r_below_in' and outer radius 'r_above_in'.
------------------------------------------------------------------------------- */
/*
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
*/

AVBMCMolIn::AVBMCMolIn(System* system_in, Box* box_in,
                       std::vector<Particle> particles_in,
                       const double r_max_inner_in,
                       const double r_below_in, const double r_above_in)
    : Moves(system_in)
{
    box = box_in;
    r_below = r_below_in;
    r_above = r_above_in;
    r_max_inner = r_max_inner_in;
    r_abovesq = r_above * r_above;
    r_belowsq = r_below * r_below;
    particles = particles_in;
    natom = particles.size();
    v_in = 1.; // 4 * pi * std::pow(r_above, 3)/3;
    label = "AVBMCMolIn";

    // ensure that first particle is located at origin
    for (Particle &particle : particles) {
        particle.r -= particles[0].r;
        particle.type = system->forcefield->label2type.at(particle.label);
    }
}


/* ----------------------------------------------------------------------------
   Insert molecule into the bonded region of an equivalent molecule
------------------------------------------------------------------------------- */

void AVBMCMolIn::perform_move()
{
    // pick molecule type to be inserted
    //MoleculeTypes* molecule_types = system->molecule_types;
    //int mol_idx = rng->choice(molecule_types->molecule_probs);
    //std::vector<std::valarray<double> > positions = molecule_types->default_mols[mol_idx];
    //natom = molecule_types->molecule_elements[mol_idx].size();

    // rotate molecule arbitrary
    //positions = rotate_molecule(atom_positions);

    // detect target molecule
    bool detected = false;
    std::vector<int> target_molecule = detect_molecule(box->particles, particles, detected, r_max_inner);
    //std::vector<int> target_molecule = system->molecule_types->detect_molecule(box->particles, mol_idx, detected);
    if (!detected) {
        reject_move = true;
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
        particles = rotate_molecule(particles);
        for (unsigned int j=0; j < natom; j++) {
            particles[j].r += box->particles[i].r + dr;
            box->particles.push_back(particles[j]);
            box->npar ++;
        }

        // compute energy difference
        du = 0.;
        for (int j=0; j < natom; j++) {
            du += system->forcefield->comp_energy_par(box->particles, box->npar - j - 1);
        }
        box->poteng += du;
    }
}


/* ----------------------------------------------------------------------------
   Get acceptance probability of move
------------------------------------------------------------------------------- */

double AVBMCMolIn::accept(double temp, double chempot)
{
    if (reject_move) {
        return 0.;
    }
    else {
        bool accept_boundary = box->boundary->correct_position();
        double dw = system->sampler->w(box->npar) - system->sampler->w(box->npar - natom);
        double prefactor = (v_in * box->npar) / ((nmolavg + 1) * (box->npar + natom)); 
        return prefactor * accept_boundary * std::exp(-(du-chempot+dw)/temp);
    }
}


/* ----------------------------------------------------------------------------
   Set back to old state before move
------------------------------------------------------------------------------- */

void AVBMCMolIn::reset()
{
    if (!reject_move) {
        box->npar -= natom;
        box->poteng -= du;
        box->particles.erase(box->particles.begin() + box->npar);
    }
}


/* ----------------------------------------------------------------------------
   Update number of time this system size has occured if
   move was accepted
------------------------------------------------------------------------------- */

void AVBMCMolIn::update_nsystemsize()
{
    if (box->npar + 1 > box->nsystemsize.size()) {
        box->nsystemsize.resize(box->npar + 1);
    }
    box->nsystemsize[box->npar] ++;
}


/* ----------------------------------------------------------------------------
   Represent move in a clean way
------------------------------------------------------------------------------- */

std::string AVBMCMolIn::repr()
{
    std::string move_info;
    move_info += "AVBMC molecule insertion move\n";
    move_info += "    Radius of outer sphere: " + std::to_string(r_above) + "\n";
    move_info += "    Radius of inner sphere: " + std::to_string(r_below) + "\n";
    return move_info;
}
