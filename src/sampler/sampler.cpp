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
    // sample nmoves moves
    for (int i=0; i<nmoves; i++){
        // pick move type
        for (double prob : system->moves_prob) {
            std::cout << "prob sample " << prob << std::endl;
        }
        std::cout << "sample1" << std::endl;
        Moves* move = system->moves[rng->choice(system->moves_prob)];
        std::cout << "sample2" << std::endl;
        move->ndrawn ++;
        std::cout << "sample3" << std::endl;
        move->perform_move();
        std::cout << "sample4" << std::endl;
        if(move->accept(system->temp, system->chempot) > rng->next_double()) {
            move->naccept ++;
        }
        else {
            move->reset();
        }
        std::cout << "sample5" << std::endl;
        move->update_nsystemsize();
        std::cout << "sample6" << std::endl;
    }
}
