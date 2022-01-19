/* ----------------------------------------------------------------------------
------------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------
   Energy-biased Aggregation-volume-biased Monte Carlo (EB-AVBMC) insertion
   and deletion moves. Here, 50% of the attemped moves are insertion moves,
   and 50% of the moves are deletion moves to maintain detailed balance. The
   EB-AVBMC type of moves where first proposed by Loeffler, Sepehri and Chen
   (2015).
------------------------------------------------------------------------------- */

#include <iostream>

#include "ebavbmc.h"
#include "../box.h"
#include "../system.h"
#include "../rng/rng.h"


/* ----------------------------------------------------------------------------
   AVBMC constructor
------------------------------------------------------------------------------- */

EBAVBMC::EBAVBMC(System* system_in, Box* box_in, const std::string label_in,
             const double r_below_in, const double r_above_in)
    : Moves(system_in), EBAVBMCIn(system_in, box_in, label_in, r_below_in, r_above_in),
      EBAVBMCOut(system_in, box_in, label_in, r_above_in)
{
    box = box_in;
    r_below = r_below_in;
    r_above = r_above_in;
    particle_label = label_in;
    label = "EB-AVBMC";
}


/* ----------------------------------------------------------------------------
   Pick in or out moves with the same probability
------------------------------------------------------------------------------- */

void EBAVBMC::perform_move()
{
    if(box->npar < 2){
        move_in = true;
    }
    else{
        move_in = rng->next_int(2);
    }
    if(move_in){
        EBAVBMCIn::perform_move();
    }
    else{
        EBAVBMCOut::perform_move();
    }
}


/* ----------------------------------------------------------------------------
   Get acceptance probability
------------------------------------------------------------------------------- */

double EBAVBMC::accept(const double temp, const double chempot)
{
    if(move_in){
        return EBAVBMCIn::accept(temp, chempot);
    }
    else{
        return EBAVBMCOut::accept(temp, chempot);
    }
}


/* ----------------------------------------------------------------------------
   Set back to old state before move as performed
------------------------------------------------------------------------------- */

void EBAVBMC::reset()
{
    if(move_in){
        EBAVBMCIn::reset();
    }
    else{
        EBAVBMCOut::reset();
    }
}


/* ----------------------------------------------------------------------------
   Update number of time this system size has occured if
   move was accepted
------------------------------------------------------------------------------- */

void EBAVBMC::update_nsystemsize()
{
    if (move_in) {
        EBAVBMCIn::update_nsystemsize();
    }
    else {
        EBAVBMCOut::update_nsystemsize();
    }
}


/* ----------------------------------------------------------------------------
   Represent move in a clean way
------------------------------------------------------------------------------- */

std::string EBAVBMC::repr()
{
    std::string move_info;
    move_info += "EB-AVBMC particle moves\n";
    move_info += "    Radius of outer sphere: " + std::to_string(r_above) + "\n";
    move_info += "    Radius of inner sphere: " + std::to_string(r_below) + "\n";
    move_info += "    Label of inserted atom: " + particle_label + "\n";
    return move_info;
}
