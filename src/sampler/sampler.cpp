#include <iostream>
#include <vector>

#include "sampler.h"
#include "../system.h"
#include "../rng/rng.h"
#include "../moves/moves.h"


/* ------------------------------------------------------
   Sampler constructor, setting box and random number 
   generator
--------------------------------------------------------- */

Sampler::Sampler(System* system_in)
{
    system = system_in;
    rng = system->rng;
}


/* -------------------------------------------------------
   Sample 'nmoves' random moves
---------------------------------------------------------- */

void Sampler::sample(int nmoves)
{
    for (int i=0; i<nmoves; i++){
        Moves* move = system->moves[rng->choice(system->moves_prob)];
        move->ndrawn ++;
        move->perform_move();
        if(move->accept(system->temp, system->chempot) > rng->next_double()) {
            move->naccept ++;
        }
        else {
            move->reset();
        }
        move->update_nsystemsize();
    }
}
