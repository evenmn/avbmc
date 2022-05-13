/* ----------------------------------------------------------------------------
   AVBMC swap move from one box to another. 
------------------------------------------------------------------------------- */

#include <string>

#include "avbmcin.h"
#include "avbmcout.h"
#include "avbmcswapright.h"
#include "../system.h"
#include "../box.h"


/* ----------------------------------------------------------------------------
   AVBMC Swap Right class constructor. Moving a particle from box 1 to box 2.
   Number of particles is conserved in the system, but not in the two boxes.
------------------------------------------------------------------------------- */

AVBMCSwapRight::AVBMCSwapRight(System *system_in, Box *box1_in, Box *box2_in,
    const std::string &particle_in, double r_below_in, double r_above_in,
    bool energy_bias_in)
    : Moves(system_in),
      AVBMCIn(system_in, box2_in, particle_in, r_below_in, r_above_in, energy_bias_in),
      AVBMCOut(system_in, box1_in, particle_in, r_above_in, energy_bias_in)
{
    box1 = box1_in;
    box2 = box2_in;
    energy_bias = energy_bias_in;
    r_below = r_below_in;
    r_above = r_above_in;
    label = "AVBMCSwapRight";
}


/* ----------------------------------------------------------------------------
   Perform swap move
------------------------------------------------------------------------------- */

void AVBMCSwapRight::perform_move()
{
    AVBMCOut::perform_move();  // remove particle from box 1
    AVBMCIn::perform_move();   // insert particle into box 2
}


/* ----------------------------------------------------------------------------
   Accept move
------------------------------------------------------------------------------- */

double AVBMCSwapRight::accept(double temp, double chempot)
{
    double acceptout = AVBMCOut::accept(temp, chempot);
    double acceptin = AVBMCIn::accept(temp, chempot);
    return acceptout * acceptin;
}


/* ----------------------------------------------------------------------------
   Reset if move was rejected
------------------------------------------------------------------------------- */

void AVBMCSwapRight::reset()
{
    AVBMCOut::reset();
    AVBMCIn::reset();
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
