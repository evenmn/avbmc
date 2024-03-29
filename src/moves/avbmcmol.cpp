/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.

  Author(s): Even M. Nordhagen
  Email(s): evenmn@mn.uio.no
  Date: 2022-06-03 (last changed 2022-08-28)
---------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------
  AVBMC molecule insertion and deletion moves. The AVBMC moves were first
  proposed by Chen and Siepmann, "A Novel Monte Carlo Algorithm for Simulating
  Strongly Associating Fluids:  Applications to Water, Hydrogen Fluoride, and
  Acetic Acid" (2000). Non-bonded molecule moves have not yet been reported
  in the literature
---------------------------------------------------------------------------- */

#include <iostream>
#include <string>
#include <vector>

#include "avbmcmol.h"
#include "avbmcmolin.h"
#include "avbmcmolout.h"

#include "../box.h"
#include "../system.h"
#include "../particle.h"
#include "../rng/rng.h"


/* ----------------------------------------------------------------------------
   AVBMC Molecule class constructor, which performs 50% insertation moves and
   50% deletion moves. 
---------------------------------------------------------------------------- */

AVBMCMol::AVBMCMol(System* system_in, Box* box_in, std::vector<Particle> molecule_in,
    const double r_below_in, const double r_above_in, const double r_inner_in,
    bool energy_bias_in, bool target_mol_in)
    : Moves(system_in),
      AVBMCMolIn(system_in, box_in, molecule_in, r_below_in, r_above_in, r_inner_in),
      AVBMCMolOut(system_in, box_in, molecule_in, r_above_in, r_inner_in)
{
    box = box_in;
    energy_bias = energy_bias_in;
    target_mol = target_mol_in;
    r_below = r_below_in;
    r_above = r_above_in;
    r_inner = r_inner_in;
    label = "AVBMCMol";
}


/* ----------------------------------------------------------------------------
   Pick in or out moves with the same probability
---------------------------------------------------------------------------- */

void AVBMCMol::perform_move()
{
    if (box->npar < 2) {
        move_in = true;
    }
    else {
        move_in = rng->next_int(2);
    }
    if (move_in) {
        AVBMCMolIn::perform_move();
    }
    else {
        AVBMCMolOut::perform_move();
    }
}


/* ----------------------------------------------------------------------------
   Get acceptance probability
---------------------------------------------------------------------------- */

double AVBMCMol::accept(const double beta, const double chempot)
{
    if (move_in) {
        return AVBMCMolIn::accept(beta, chempot);
    }
    else {
        return AVBMCMolOut::accept(beta, chempot);
    }
}


/* ----------------------------------------------------------------------------
   Set back to old state before move as performed
---------------------------------------------------------------------------- */

void AVBMCMol::reset()
{
    if (move_in) {
        AVBMCMolIn::reset();
    }
    else {
        AVBMCMolOut::reset();
    }
}


/* ----------------------------------------------------------------------------
   Update number of time this system size has occured if move was accepted
---------------------------------------------------------------------------- */

void AVBMCMol::update_size_histogram()
{
    box->update_size_histogram();
}


/* ----------------------------------------------------------------------------
   Represent move in a clean way
---------------------------------------------------------------------------- */

std::string AVBMCMol::repr()
{
    std::string move_info;
    move_info += "AVBMC molecule moves\n";
    move_info += "    Radius of outer sphere: " + std::to_string(r_above) + "\n";
    move_info += "    Radius of inner sphere: " + std::to_string(r_below) + "\n";
    return move_info;
}
