/* ----------------------------------------------------------------------------
------------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------
   Aggregation-volume-biased Monte Carlo (AVBMC) insertion and deletion moves.
   Here, 50% of the attemped moves are insertion moves, and 50% of the moves
   are deletion moves to maintain detailed balance. The AVBMC type of moves
   where first proposed by Chen (2000).
------------------------------------------------------------------------------- */

#include <iostream>

#include "avbmc.h"
#include "../box.h"
#include "../system.h"
#include "../rng/rng.h"


/* ----------------------------------------------------------------------------
   AVBMC constructor
------------------------------------------------------------------------------- */

AVBMC::AVBMC(System* system_in, Box* box_in, const std::string label_in,
             const double r_below_in, const double r_above_in)
    : Moves(system_in), AVBMCIn(system_in, box_in, label_in, r_below_in, r_above_in),
      AVBMCOut(system_in, box_in, label_in, r_above_in)
{
    box = box_in;
    r_below = r_below_in;
    r_above = r_above_in;
    particle_label = label_in;
    label = "AVBMC";
}


/* ----------------------------------------------------------------------------
   Pick in or out moves with the same probability
------------------------------------------------------------------------------- */

void AVBMC::perform_move()
{
    if(box->npar < 2){
        move_in = true;
    }
    else{
        move_in = rng->next_int(2);
    }
    if(move_in){
        AVBMCIn::perform_move();
    }
    else{
        AVBMCOut::perform_move();
    }
}


/* ----------------------------------------------------------------------------
   Get acceptance probability
------------------------------------------------------------------------------- */

double AVBMC::accept(const double temp, const double chempot)
{
    if(move_in){
        return AVBMCIn::accept(temp, chempot);
    }
    else{
        return AVBMCOut::accept(temp, chempot);
    }
}


/* ----------------------------------------------------------------------------
   Set back to old state before move as performed
------------------------------------------------------------------------------- */

void AVBMC::reset()
{
    if(move_in){
        AVBMCIn::reset();
    }
    else{
        AVBMCOut::reset();
    }
}


/* ----------------------------------------------------------------------------
   Update number of time this system size has occured if
   move was accepted
------------------------------------------------------------------------------- */

void AVBMC::update_nsystemsize()
{
    if (move_in) {
        AVBMCIn::update_nsystemsize();
    }
    else {
        AVBMCOut::update_nsystemsize();
    }
}


/* ----------------------------------------------------------------------------
   Represent move in a clean way
------------------------------------------------------------------------------- */

std::string AVBMC::repr()
{
    std::string move_info;
    move_info += "AVBMC particle moves\n";
    move_info += "    Radius of outer sphere: " + std::to_string(r_above) + "\n";
    move_info += "    Radius of inner sphere: " + std::to_string(r_below) + "\n";
    move_info += "    Label of inserted atom: " + particle_label + "\n";
    return move_info;
}
