#include <iostream>
#include <cmath>
#include <valarray>
#include <vector>

#include "avbmcinmol.h"
#include "../box.h"
#include "../particle.h"
#include "../molecule.h"


/* -----------------------------------------------------
   Aggregate-volume-biased in move. This is a fictitious 
   inter-box move, and works in the grand canonical 
   essemble only. A defined molecule is inserted
-------------------------------------------------------- */


/* -----------------------------------------------------
   AVBMC insertation mol move constructor. Molecule is
   inserted in the bonded region of an equivalent
   molecule, defined by inner radius 'r_below_in' and
   outer radius 'r_above_in'. 
-------------------------------------------------------- */

AVBMCInMol::AVBMCInMol(Box* box_in, const double r_below_in, const double r_above_in)
    : Moves(box_in)
{
    r_below = r_below_in;
    r_above = r_above_in;
    r_abovesq = r_above * r_above;
    r_belowsq = r_below * r_below;
    v_in = 1.; // 4 * pi * std::pow(r_above, 3)/3; // can be set to 1 according to Henrik
}


/* ------------------------------------------------------------
   Insert molecule into the bonded region of an equivalent
   molecule
--------------------------------------------------------------- */

void AVBMCInMol::perform_move()
{
    // pick molecule type to be inserted
    int mol_idx = box->rng->choice(box->molecule_types->molecule_probs);
    std::vector<std::valarray<double> > positions = box->molecule_types->default_mols[mol_idx];
    natom = box->molecule_types->molecule_elements[mol_idx].size();
    //int com_type = box->molecule_types->molecule_types[mol_idx][0];

    // rotate molecule arbitrary
    positions = rotate_molecule(positions);

    // search for target molecule
    bool constructed = false;
    Molecule* molecule_target = box->molecule_types->construct_molecule(mol_idx, box->particles, constructed);
    if (!constructed) {
        reject_move = true;
    }
    else {
        reject_move = false;
        int i = molecule_target->atoms_idx[0];  // center atom
        std::vector<int> neigh_listi = box->forcefield->build_neigh_list(i, r_abovesq);
        int n_in = neigh_listi.size();
        nmolavg = (double) n_in / natom;

    /*
    // search for COM particle (first particle)
    int type, i;
    int count = 0;
    while(type != com_type && count < box->npar){
        i = box->rng->next_int(box->npar);
        type = box->particles[i]->type;
        count ++;
    }

    if (type == com_type){
        reject_move = false;
    */
        // compute norm
        /*
        auto norm = [] (std::valarray<double> x) -> double { 
            double sqrd_sum = 0.;
            for(double x_ : x){
                sqrd_sum += x_ * x_;
            }
            return sqrd_sum;
        };
        */

        // shift out molecule relative to target molecule
        std::valarray<double> dr(box->ndim);
        double normsq = std::pow(dr, 2).sum(); //norm(dr);
        while(normsq > r_abovesq || normsq < r_belowsq){
            for(double &d : dr){
                d = r_above * (2 * box->rng->next_double() - 1);
            }
            normsq = std::pow(dr, 2).sum();
        }

        // construct new particles
        for(int j=0; j < natom; j++){
            positions[j] += box->particles[i]->r + dr;
            std::string element = box->molecule_types->molecule_elements[mol_idx][j];
            Particle* particle = new Particle(element, positions[j]);
            for(int k=0; k<box->ntype; k++){
                if (box->unique_labels[k] == element){
                    particle->type = k;
                }
            }
            box->particles.push_back(particle);
            box->npar ++;
        }

        // compute energy difference
        du = 0.;
        for(int j=0; j < natom; j++){
            du += box->forcefield->comp_energy_par(box->particles, box->npar - j - 1);
        }
    }
}


/* ---------------------------------------------------------
   Get acceptance probability of move
------------------------------------------------------------ */

double AVBMCInMol::accept(double temp, double chempot)
{
    if (reject_move) {
        return 0.;
    }
    else {
        double dw = box->sampler->w(box->npar) - box->sampler->w(box->npar - natom);
        return (v_in * box->npar) / ((nmolavg + 1) * (box->npar + natom)) * std::exp(-(du-chempot+dw)/temp);
    }
}


/* ----------------------------------------------------------
   Set back to old state before move
------------------------------------------------------------- */

void AVBMCInMol::reset()
{
    if (!reject_move) {
        box->npar -= natom;
        box->poteng -= du;
        box->particles.erase(box->particles.begin() + box->npar);
    }
}


/* -----------------------------------------------------
   Represent move in a clean way
-------------------------------------------------------- */

std::string AVBMCInMol::repr()
{
    std::string move_info;
    move_info += "    Label: AVBMC insertion molecule move\n";
    move_info += "    Radius of outer sphere: " + std::to_string(r_above) + "\n";
    move_info += "    Radius of inner sphere: " + std::to_string(r_below) + "\n";
    return move_info;
}
