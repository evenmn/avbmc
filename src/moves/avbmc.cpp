#include <iostream>

#include "avbmc.h"
#include "../box.h"


/* -----------------------------------------------------------
   AVBMC constructor, which performs 50% insertation moves and
   50% deletion moves. 
-------------------------------------------------------------- */

AVBMC::AVBMC(Box* box_in, const double r_below_in, const double r_above_in)
    : AVBMCIn(box_in, r_below_in, r_above_in), AVBMCOut(box_in, r_above_in), Moves(box_in) {}


/* -----------------------------------------------------------
   Pick in or out moves with the same probability
-------------------------------------------------------------- */

void AVBMC::perform_move()
{
    if(AVBMCIn::box->npar < 2){
        move_in = true;
    }
    else{
        move_in = AVBMCIn::box->rng->next_int(2);
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
