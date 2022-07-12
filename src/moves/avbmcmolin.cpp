/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.

  Author(s): Even M. Nordhagen
  Email(s): evenmn@mn.uio.no
  Date: 2022-06-03 (last changed 2022-06-03)
---------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------
  AVBMC molecule insertion moves. The AVBMC moves were first proposed by Chen
  and Siepmann, "A Novel Monte Carlo Algorithm for Simulating Strongly
  Associating Fluids:  Applications to Water, Hydrogen Fluoride, and Acetic
  Acid" (2000). Non-bonded molecule moves have not yet been reported in the
  literature
---------------------------------------------------------------------------- */

#include <iostream>
#include <cmath>
#include <valarray>
#include <vector>
#include <cassert>

#include "avbmcmolin.h"
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
   AVBMC molecule insertation move constructor. Molecule is inserted in the
   bonded region of an equivalent molecule, defined by inner radius
   'r_below_in' and outer radius 'r_above_in'.
---------------------------------------------------------------------------- */

AVBMCMolIn::AVBMCMolIn(System* system_in, Box* box_in,
    std::vector<Particle> molecule_in, const double r_below_in,
    const double r_above_in, const double r_inner_in, const bool energy_bias_in,
    const bool target_mol_in)
    : Moves(system_in)
{
    box = box_in;
    r_inner = r_inner_in;
    r_below = r_below_in;
    r_above = r_above_in;
    r_abovesq = r_above * r_above;
    r_belowsq = r_below * r_below;

    energy_bias = energy_bias_in;
    target_mol = target_mol_in;
    particles = molecule_in;
    //molecule_template = molecule_in;
    natom = particles.size();
    natom_inv = 1. / natom;
    v_in = 1.; // 4 * pi * std::pow(r_above, 3)/3;
    cum_time = 0.;
    naccept = ndrawn = 0;
    label = "AVBMCMolIn";

    /*
    neigh_id_above = box->distance_manager->add_cutoff(
        r_above, particles[0].label, particles[0].label, false);
    neigh_id_inner = box->distance_manager->add_cutoff(
        r_inner, particles[0].label, particles[1].label, false);
    */
    // ensure that first particle is located at origin
    for (Particle &particle : particles) {
        particle.r -= particles[0].r;
        particle.type = box->forcefield->label2type.at(particle.label);
    }
}


/* ----------------------------------------------------------------------------
   Detect target molecule (or atom if target_mol=false). Make maximum 'npar' 
   attempts of finding target molecule. For performance reasons, this is done
   before insertion molecule is created.
---------------------------------------------------------------------------- */
/*
unsigned int AVBMCMolIn::detect_target_molecule(bool &detected)
{
    unsigned int count, i;
    std::vector<unsigned int> target_molecule;
    std::vector<std::vector<unsigned int> > neigh_list_inner;

    count = i = 0;
    detected = false;
    while (count < box->npar && !detected) {
        count ++;
        neigh_list_inner = box->distance_manager->neigh_lists[neigh_id_inner];
        target_molecule = detect_molecule(neigh_list_inner, particles, detected);
        i = target_molecule[0];
    }
    return i;
}
*/

/* ----------------------------------------------------------------------------
   Create insertion molecule.
---------------------------------------------------------------------------- */

std::vector<Particle> AVBMCMolIn::create_molecule()
{
    std::vector<Particle> molecule;

    for (Particle particle : particles) {
        molecule.push_back(Particle(particle)); // utilizing copy constructor
    }
    molecule = rotate_molecule(molecule);
    return molecule;
}


/* ----------------------------------------------------------------------------
   Insert molecule if target can be detected
---------------------------------------------------------------------------- */

void AVBMCMolIn::insert(std::vector<Particle> molecule)
{
    unsigned int i, j;
    std::valarray<double> dr;

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
    i = box->typeidx[particles[0].type][rng->next_int(box->npartype[particles[0].type])];
    std::vector<unsigned int> neigh_listi = box->distance_manager->build_neigh_list(i, r_abovesq);
    nmolavg = neigh_listi.size() * natom_inv;

    // add molecule
    dr = box->particles[i].r + insertion_position(false);
    for (Particle &particle : molecule) {
        particle.r += dr;
        box->add_particle(particle);
    }

    // compute energy difference
    du = 0.;
    for (j=0; j < natom; j++) {
        du += box->forcefield->comp_energy_par_force0(box->npar - j - 1);
    }
    box->poteng += du;
}


/* ----------------------------------------------------------------------------
   Insert molecule into the bonded region of an equivalent molecule
---------------------------------------------------------------------------- */
        
void AVBMCMolIn::perform_move()
{
    insert(create_molecule());
}


/* ----------------------------------------------------------------------------
   Get acceptance probability of move
---------------------------------------------------------------------------- */

double AVBMCMolIn::accept(double temp, double chempot)
{
    for (Constraint* constraint : box->constraints) {
        if (!constraint->verify()) return 0.;
    }
    double dw = system->sampler->w(box->npar) - system->sampler->w(box->npar - natom);
    double prefactor = (v_in * box->npar) / ((nmolavg + 1) * (box->npar + natom)); 
    double accept_prob = prefactor * std::exp(-(du-chempot+dw)/temp);
    return accept_prob;
}


/* ----------------------------------------------------------------------------
   Set back to old state before move
---------------------------------------------------------------------------- */

void AVBMCMolIn::reset()
{
    for (unsigned int i=natom; i--;) {
        box->rm_particle(box->npar-1);
    }
    box->poteng -= du;
}


/* ----------------------------------------------------------------------------
   Update number of time this system size has occured if move was accepted
---------------------------------------------------------------------------- */

void AVBMCMolIn::update_size_histogram()
{
    box->update_size_histogram();
}


/* ----------------------------------------------------------------------------
   Represent move in a clean way
---------------------------------------------------------------------------- */

std::string AVBMCMolIn::repr()
{
    std::string move_info;
    move_info += "AVBMC molecule insertion move\n";
    move_info += "    Radius of outer sphere: " + std::to_string(r_above) + "\n";
    move_info += "    Radius of inner sphere: " + std::to_string(r_below) + "\n";
    move_info += "    Energy bias:            " + std::to_string(energy_bias) + "\n";
    move_info += "    Search target molecule: " + std::to_string(target_mol) + "\n";
    return move_info;
}
