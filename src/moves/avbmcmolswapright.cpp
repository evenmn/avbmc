/* ----------------------------------------------------------------------------
   AVBMC swap move from one box to another. 
---------------------------------------------------------------------------- */

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
   boxes separately.
---------------------------------------------------------------------------- */

AVBMCMolSwapRight::AVBMCMolSwapRight(System *system_in, Box *box1_in,
    Box *box2_in, std::vector<Particle> molecule_in, double r_below_in,
    double r_above_in, double r_inner_in, bool energy_bias_in,
    bool target_mol_in)
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
---------------------------------------------------------------------------- */

void AVBMCMolSwapRight::perform_move()
{
    AVBMCMolOut::perform_move();  // remove molecule from box 1

    // do not attempt inserting molecule if molecule was not removed
    if (!AVBMCMolOut::detected_out) return;

    
    AVBMCMolIn::insert(AVBMCMolOut::molecule_out);
}


/* ----------------------------------------------------------------------------
   Acceptance probability
---------------------------------------------------------------------------- */

double AVBMCMolSwapRight::accept(double temp, double chempot)
{
    double acceptin = AVBMCMolIn::accept(temp, chempot);
    double acceptout = AVBMCMolOut::accept(temp, chempot);
    return std::min(acceptin, acceptout);
}


/* ----------------------------------------------------------------------------
   Reset if move was rejected
---------------------------------------------------------------------------- */

void AVBMCMolSwapRight::reset()
{
    if (AVBMCMolOut::detected_out) {
        AVBMCMolIn::reset();  // this order is important: always reset
        AVBMCMolOut::reset(); // AVBMCMolIn before AVBMCMolOut
    }
}


/* ----------------------------------------------------------------------------
   Update number of time this system size has occured if move was accepted
---------------------------------------------------------------------------- */

void AVBMCMolSwapRight::update_size_histogram()
{
    box1->update_size_histogram();
    box2->update_size_histogram();
}


/* ----------------------------------------------------------------------------
   Represent move
---------------------------------------------------------------------------- */

std::string AVBMCMolSwapRight::repr()
{
    return "";
}
