#include <iostream>

#include "avbmc.h"
#include "../box.h"
#include "../system.h"
#include "../rng/rng.h"


/* -----------------------------------------------------------
   AVBMC constructor, which performs 50% insertation moves and
   50% deletion moves. 
-------------------------------------------------------------- */

AVBMC::AVBMC(System* system_in, Box* box_in, const double r_below_in, const double r_above_in)
    : AVBMCIn(system_in, box_in, r_below_in, r_above_in),
      AVBMCOut(system_in, box_in, r_above_in), Moves(system_in) 
{
    box = box_in;
    boxes.push_back(box);
    r_below = r_below_in;
    r_above = r_above_in;
    label = "AVBMC";
}


/* -----------------------------------------------------------
   Pick in or out moves with the same probability
-------------------------------------------------------------- */

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


/* -----------------------------------------------------------
   Get acceptance probability
-------------------------------------------------------------- */

double AVBMC::accept(const double temp, const double chempot)
{
    if(move_in){
        return AVBMCIn::accept(temp, chempot);
    }
    else{
        return AVBMCOut::accept(temp, chempot);
    }
}


/* -----------------------------------------------------------
   Set back to old state before move as performed
-------------------------------------------------------------- */

void AVBMC::reset()
{
    if(move_in){
        AVBMCIn::reset();
    }
    else{
        AVBMCOut::reset();
    }
}


/* -----------------------------------------------------
   Represent move in a clean way
-------------------------------------------------------- */

std::string AVBMC::repr()
{
    std::string move_info;
    move_info += "AVBMC particle moves\n";
    move_info += "    Radius of outer sphere: " + std::to_string(r_above) + "\n";
    move_info += "    Radius of inner sphere: " + std::to_string(r_below) + "\n";
    return move_info;
}
