#include <iostream>
#include <string>

#include "avbmcmol.h"
#include "avbmcmolin.h"
#include "avbmcmolout.h"
#include "../box.h"
#include "../system.h"
#include "../rng/rng.h"


/* -----------------------------------------------------------
   AVBMC constructor, which performs 50% insertation moves and
   50% deletion moves. 
-------------------------------------------------------------- */

AVBMCMol::AVBMCMol(System* system_in, Box* box_in, std::vector<class Particle> particles_in, double r_inner_max_in, const double r_below_in, const double r_above_in)
    : Moves(system_in), AVBMCMolIn(system_in, box_in, particles_in, r_below_in, r_above_in, r_inner_max_in), AVBMCMolOut(system_in, box_in, r_above_in)
{
    box = box_in;
    r_below = r_below_in;
    r_above = r_above_in;
    label = "AVBMCMol";
}


/* -----------------------------------------------------------
   Pick in or out moves with the same probability
-------------------------------------------------------------- */

void AVBMCMol::perform_move()
{
    if(box->npar < 2){
        move_in = true;
    }
    else{
        move_in = rng->next_int(2);
    }
    if(move_in){
        AVBMCMolIn::perform_move();
    }
    else{
        AVBMCMolOut::perform_move();
    }
}


/* -----------------------------------------------------------
   Get acceptance probability
-------------------------------------------------------------- */

double AVBMCMol::accept(const double temp, const double chempot)
{
    if(move_in){
        return AVBMCMolIn::accept(temp, chempot);
    }
    else{
        return AVBMCMolOut::accept(temp, chempot);
    }
}


/* -----------------------------------------------------------
   Set back to old state before move as performed
-------------------------------------------------------------- */

void AVBMCMol::reset()
{
    if(move_in){
        AVBMCMolIn::reset();
    }
    else{
        AVBMCMolOut::reset();
    }
}


/* ----------------------------------------------------------
   Update number of time this system size has occured if
   move was accepted
------------------------------------------------------------- */

void AVBMCMol::update_nsystemsize()
{
    if (move_in) {
        AVBMCMolIn::update_nsystemsize();
    }
    else {
        AVBMCMolOut::update_nsystemsize();
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
