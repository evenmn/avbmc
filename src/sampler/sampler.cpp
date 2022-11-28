/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.

  Author(s): Even M. Nordhagen
  Email(s): evenmn@mn.uio.no
  Date: 2022-06-03 (last changed 2022-06-03)
---------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------
  File hosting the sampler base class. Used to perform a Monte Carlo cycle
---------------------------------------------------------------------------- */

#include <iostream>
#include <vector>
#include <chrono>

#include "sampler.h"
#include "../system.h"
#include "../rng/rng.h"
#include "../moves/moves.h"


/* ----------------------------------------------------------------------------
   Sampler constructor, setting box and random number generator.
---------------------------------------------------------------------------- */

Sampler::Sampler(System* system_in)
{
    system = system_in;
    rng = system->rng;
    bool update_histogram = false;
}


/* ----------------------------------------------------------------------------
   Sample 'nmoves' random moves.
---------------------------------------------------------------------------- */

void Sampler::sample(int nmoves)
{
    //typedef std::chrono::high_resolution_clock Time;
    //typedef std::chrono::duration<double> fsec;
    //auto t0 = Time::now();
    for (int i=0; i<nmoves; i++){
        int j = rng->choice(system->moves_prob);
        Moves* move = system->moves[j];
        move->ndrawn ++;
        move->perform_move();
        if (move->accept(system->beta, system->chempot) > rng->next_double()) {
            move->naccept ++;
        }
        else {
            move->reset();
        }
        //if (update_histogram) move->update_size_histogram();
        //auto t1 = Time::now();
        //fsec fs = t1 - t0;
        //move->cum_time += fs.count();
        //t0 = t1;
    }
}
