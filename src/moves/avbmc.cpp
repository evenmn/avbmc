#include "avbmc.h"
#include "../box.h"


AVBMC::AVBMC(class Box* box_in, const double r_below_in, const double r_above_in)
    : AVBMCIn(box_in, r_below_in, r_above_in), AVBMCOut(box_in, r_above_in), Moves(box_in) {}

void AVBMC::perform_move()
{
    /* Pick in or out with the same probability
     */

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

double AVBMC::accept()
{
    /* 
     */
    if(move_in){
        return AVBMCIn::accept();
    }
    else{
        return AVBMCOut::accept();
    }
}

void AVBMC::reset()
{
    if(move_in){
        AVBMCIn::reset();
    }
    else{
        AVBMCOut::reset();
    }
}
