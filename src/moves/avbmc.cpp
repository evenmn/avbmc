#include "avbmc.h"
#include "../box.h"


AVBMC::AVBMC(class Box* box_in, const double r_below_in, const double r_above_in)
    : AVBMCIn(box_in, r_below_in, r_above_in), AVBMCOut(box_in, r_above_in), Moves(box_in) {}

void AVBMC::perform_move(const int i)
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
        AVBMCIn::perform_move(i);
    }
    else{
        AVBMCOut::perform_move(i);
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

void AVBMC::update_box(const int i)
{
    /*
     */
    if(move_in){
        AVBMCIn::update_box(i);
    }
    else{
        AVBMCOut::update_box(i);
    }
}
