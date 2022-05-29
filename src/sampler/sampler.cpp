#include <iostream>
#include <vector>
#include <chrono>

#include "sampler.h"
#include "../system.h"
#include "../rng/rng.h"
#include "../moves/moves.h"


/* ----------------------------------------------------------------------------
   Sampler constructor, setting box and random number generator.
------------------------------------------------------------------------------- */

Sampler::Sampler(System* system_in)
{
    system = system_in;
    rng = system->rng;
}


/* ----------------------------------------------------------------------------
   Sample 'nmoves' random moves.
------------------------------------------------------------------------------- */

void Sampler::sample(int nmoves)
{
    int i, j;
    typedef std::chrono::high_resolution_clock Time;
    typedef std::chrono::duration<double> fsec;
    auto t0 = Time::now();
    for (i=0; i<nmoves; i++){
        j = rng->choice(system->moves_prob);
        Moves* move = system->moves[j];
        move->ndrawn ++;
        move->perform_move();
        if(move->accept(system->temp, system->chempot) > rng->next_double()) {
            move->naccept ++;
        }
        else {
            move->reset();
        }
        move->update_size_histogram();
        auto t1 = Time::now();
        fsec fs = t1 - t0;
        move->cum_time += fs.count();
        t0 = t1;
    }
}
