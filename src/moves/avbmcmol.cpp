#include <iostream>
#include <string>

#include "avbmcmol.h"
#include "../box.h"


/* -----------------------------------------------------------
   AVBMC constructor, which performs 50% insertation moves and
   50% deletion moves. 
-------------------------------------------------------------- */

AVBMCMol::AVBMCMol(Box* box_in, const double r_below_in, const double r_above_in)
    : AVBMCInMol(box_in, r_below_in, r_above_in), AVBMCOutMol(box_in, r_above_in), Moves(box_in) 
{
    r_below = r_below_in;
    r_above = r_above_in;
}


/* -----------------------------------------------------------
   Pick in or out moves with the same probability
-------------------------------------------------------------- */

void AVBMCMol::perform_move()
{
    if(AVBMCInMol::box->npar < 2){
        move_in = true;
    }
    else{
        move_in = AVBMCInMol::box->rng->next_int(2);
    }
    if(move_in){
        AVBMCInMol::perform_move();
    }
    else{
        AVBMCOutMol::perform_move();
    }
}


/* -----------------------------------------------------------
   Get acceptance probability
-------------------------------------------------------------- */

double AVBMCMol::accept(const double temp, const double chempot)
{
    if(move_in){
        return AVBMCInMol::accept(temp, chempot);
    }
    else{
        return AVBMCOutMol::accept(temp, chempot);
    }
}


/* -----------------------------------------------------------
   Set back to old state before move as performed
-------------------------------------------------------------- */

void AVBMCMol::reset()
{
    if(move_in){
        AVBMCInMol::reset();
    }
    else{
        AVBMCOutMol::reset();
    }
}


/* -----------------------------------------------------
   Represent move in a clean way
-------------------------------------------------------- */

std::string AVBMCMol::repr()
{
    std::string move_info;
    move_info += "AVBMC molecule moves\n";
    move_info += "    Radius of outer sphere: " + std::to_string(r_above) + "\n";
    move_info += "    Radius of inner sphere: " + std::to_string(r_below) + "\n";
    return move_info;
}
