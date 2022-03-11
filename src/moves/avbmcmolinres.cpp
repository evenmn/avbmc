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

#include "avbmcmolinres.h"
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
------------------------------------------------------------------------------- */

AVBMCMolInRes::AVBMCMolInRes(System* system_in, Box* box_in,
                       std::vector<Particle> particles_in,
                       const double r_inner_in,
                       const double r_below_in, const double r_above_in,
                       const bool energy_bias_in, const bool target_mol_in)
    : Moves(system_in)
{
    box = box_in;
    r_below = r_below_in;
    r_above = r_above_in;
    r_inner = r_inner_in;
    r_abovesq = r_above * r_above;
    r_belowsq = r_below * r_below;

    energy_bias = energy_bias_in;
    target_mol = target_mol_in;
    particles = particles_in;
    natom = particles.size();
    natom_inv = 1. / natom;
    v_in = 1.; // 4 * pi * std::pow(r_above, 3)/3;
    label = "AVBMCMolIn";

    neigh_id_above = box->distance_manager->add_cutoff(r_above, particles[0].label, particles[0].label);
    neigh_id_inner = box->distance_manager->add_cutoff(r_inner, particles[0].label, particles[1].label);

    // ensure that first particle is located at origin
    for (Particle &particle : particles) {
        particle.r -= particles[0].r;
        particle.type = system->forcefield->label2type.at(particle.label);
    }
}


/* ----------------------------------------------------------------------------
   Insert molecule into the bonded region of an equivalent molecule
------------------------------------------------------------------------------- */

void AVBMCMolInRes::perform_move()
{
    bool detected_target;
    unsigned int count, i, j;
    std::vector<int> neigh_listi;
    std::vector<std::vector<int> > neigh_list_inner;

    count = 0;
    reject_move = true;
    detected_target = false;
    while (count < box->npar && reject_move) {
        count ++;
        if (target_mol) {
            neigh_list_inner = box->distance_manager->neigh_lists[neigh_id_inner];
            std::vector<int> target_molecule = detect_molecule(neigh_list_inner, particles, detected_target);
            i = target_molecule[0];
        }
        else {
            i = rng->next_int(box->npar);
            detected_target = (box->particles[i].type == particles[0].type);
        }
        if (detected_target) {
            reject_move = false;
            
            //neigh_listi = box->distance_manager->neigh_lists[neigh_id_above][i];
            neigh_listi = box->build_neigh_list(i, r_abovesq);
            nmolavg = neigh_listi.size() * natom_inv;

            // shift in-molecule relative to target molecule
            // replace with sin/cos
            std::valarray<double> dr(system->ndim);
            double normsq = norm(dr);
            while (normsq > r_abovesq || normsq < r_belowsq) {
                for (double &d : dr) {
                    d = r_above * (2 * rng->next_double() - 1);
                }
                normsq = norm(dr);
            }
            dr += box->particles[i].r;

            npartype_old = box->npartype;
            particles_old = box->particles;
            // construct new molecule
            particles = rotate_molecule(particles);
            box->distance_manager->set();
            for (j=0; j < natom; j++) {
                box->npar ++;
                std::valarray<double> new_pos = particles[j].r + dr;
                Particle particle(particles[j].label, new_pos);
                particle.type = particles[j].type;
                box->npartype[particle.type] ++;
                box->particles.push_back(particle);
                box->distance_manager->update_insert(box->npar - 1);
            }

            // compute energy difference
            du = 0.;
            for (j=0; j < natom; j++) {
                du += system->forcefield->comp_energy_par(box, box->npar - j - 1);
            }
            box->poteng += du;
        }
    }
}


/* ----------------------------------------------------------------------------
   Get acceptance probability of move
------------------------------------------------------------------------------- */

double AVBMCMolInRes::accept(double temp, double chempot)
{
    bool constr_satis = true;
    for (Constraint* constraint : box->constraints) {
        constr_satis *= constraint->verify();
    }
    double accept_prob = 0.;
    if (!reject_move && constr_satis) {
        bool accept_boundary = box->boundary->correct_position();
        double dw = system->sampler->w(box->npar) - system->sampler->w(box->npar - natom);
        double prefactor = (v_in * box->npar) / ((nmolavg + 1) * (box->npar + natom)); 
        accept_prob = prefactor * accept_boundary * std::exp(-(du-chempot+dw)/temp);
    }
    return accept_prob;
}


/* ----------------------------------------------------------------------------
   Set back to old state before move
------------------------------------------------------------------------------- */

void AVBMCMolInRes::reset()
{
    if (!reject_move) {
        box->distance_manager->reset();
        box->npar -= natom;
        box->poteng -= du;
        box->npartype = npartype_old;
        box->particles = particles_old;
    }
}


/* ----------------------------------------------------------------------------
   Update number of time this system size has occured if move was accepted
------------------------------------------------------------------------------- */

void AVBMCMolInRes::update_nsystemsize()
{
    if (box->npar + 1 > box->nsystemsize.size()) {
        box->nsystemsize.resize(box->npar + 1);
    }
    box->nsystemsize[box->npar] ++;
}


/* ----------------------------------------------------------------------------
   Represent move in a clean way
------------------------------------------------------------------------------- */

std::string AVBMCMolInRes::repr()
{
    std::string move_info;
    move_info += "AVBMC molecule insertion move\n";
    move_info += "    Radius of outer sphere: " + std::to_string(r_above) + "\n";
    move_info += "    Radius of inner sphere: " + std::to_string(r_below) + "\n";
    move_info += "    Energy bias:            " + std::to_string(energy_bias) + "\n";
    move_info += "    Search target molecule: " + std::to_string(target_mol) + "\n";
    return move_info;
}
