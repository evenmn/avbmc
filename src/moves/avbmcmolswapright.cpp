/* ----------------------------------------------------------------------------
   AVBMC swap move from one box to another. 
------------------------------------------------------------------------------- */

#include <string>

#include "avbmcmolin.h"
#include "avbmcmolout.h"
#include "avbmcmolswapright.h"
#include "../system.h"
#include "../box.h"
#include "../particle.h"


/* ----------------------------------------------------------------------------
   AVBMC Swap Molecule Right class constructor. Moving a particle from box 1 to
   box 2. Number of particles is conserved in the system, but not in the two
   boxes.
------------------------------------------------------------------------------- */

AVBMCMolSwapRight::AVBMCMolSwapRight(System *system_in, Box *box1_in,
    Box *box2_in, std::vector<Particle> molecule_in, double r_below_in,
    double r_above_in, double r_inner_in, bool energy_bias_in, bool target_mol_in)
    : Moves(system_in),
      AVBMCMolIn(system_in, box2_in, molecule_in, r_below_in, r_above_in,
        r_inner_in, energy_bias_in, target_mol_in),
      AVBMCMolOut(system_in, box1_in, molecule_in, r_above_in, r_inner_in,
        energy_bias_in, target_mol_in)
{
    box1 = box1_in;
    box2 = box2_in;
    energy_bias = energy_bias_in;
    target_mol = target_mol_in;
    r_below = r_below_in;
    r_above = r_above_in;
    r_inner = r_inner_in;
    label = "AVBMCMolSwapRight";
}


/* ----------------------------------------------------------------------------
   Perform swap move
------------------------------------------------------------------------------- */

void AVBMCMolSwapRight::perform_move()
{
    unsigned int i, j;
    std::vector<unsigned int> neigh_listi;
    std::valarray<double> dr;

    AVBMCMolOut::perform_move();  // remove molecule from box 1
    if (AVBMCMolOut::detected_out) {
        // do not attempt inserting molecule if molecule was not removed
        return;
    }
    i = AVBMCMolIn::detect_target_molecule(detected_target);

    if (AVBMCMolIn::detected_target) {
        neigh_listi = box2->build_neigh_list(i, r_abovesq);
        AVBMCMolIn::nmolavg = neigh_listi.size() * AVBMCMolIn::natom_inv;
        dr = AVBMCMolIn::insertion_position(i);
        for (Particle &particle : AVBMCMolOut::molecule_out) {
            box2->npar ++;
            particle.r += dr;
            box2->npartype[particle.type] ++;
            box2->particles.push_back(particle);
            box2->distance_manager->update_insert(box2->npar - 1);
        }
        // compute energy difference
        du = 0.;
        for (j=0; j < natom; j++) {
            du += box2->forcefield->comp_energy_par_force0(box->npar - j - 1);
        }
        box->poteng += du;
    }
}


/* ----------------------------------------------------------------------------
   Accept move
------------------------------------------------------------------------------- */

double AVBMCSwapRight::accept(double temp, double chempot)
{
    double acceptout = AVBMCMolOut::accept(temp, chempot);
    double acceptin = AVBMCMolIn::accept(temp, chempot);
    return acceptout * acceptin;
}


/* ----------------------------------------------------------------------------
   Reset if move was rejected
------------------------------------------------------------------------------- */

void AVBMCSwapRight::reset()
{
    AVBMCMolOut::reset();
    AVBMCMolIn::reset();
}


/* ----------------------------------------------------------------------------
   Update number of time this system size has occured if move was accepted
------------------------------------------------------------------------------- */

void AVBMCSwapRight::update_size_histogram()
{
    box1->update_size_histogram();
    box2->update_size_histogram();
}


/* ----------------------------------------------------------------------------
   Represent move
------------------------------------------------------------------------------- */

std::string AVBMCSwapRight::repr()
{
    return "";
}
