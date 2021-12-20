#include <iostream>
#include <vector>

#include "sampler.h"
#include "../box.h"
#include "../system.h"
#include "../rng/rng.h"
#include "../moves/moves.h"
#include "../boundary/boundary.h"


/* ------------------------------------------------------
   Sampler constructor, setting box and random number 
   generator
--------------------------------------------------------- */

Sampler::Sampler(System* system_in)
{
    system = system_in;
    rng = system->rng;
    acceptance_ratio = 0.;
}


/* ------------------------------------------------------
   Given a vector of moves and a vector of corresponding 
   probabilities, returning move
--------------------------------------------------------- */
/*
Moves* Sampler::propose_move(std::vector<Moves*> moves, std::vector<double> moves_prob)
{
    move_idx = rng->choice(moves_prob);
    return moves[move_idx];
}
*/

/* ------------------------------------------------------
   Decide if move should be accepted or rejected. This
   is done in two steps:
     1. Check if boundary condition is satisfied
     2. Check if acceptance criterion is satisfied
--------------------------------------------------------- */

bool Sampler::accept_move(Moves* move, double temp, double chempot)
{
    // check boundary condition
    bool accept_boundary = true;
    for(Box* box : move->boxes){
        accept_boundary *= box->boundary->correct_position();
    }

    // check move acceptance probability
    bool accept_sampler = move->accept(temp, chempot) > rng->next_double();

    // utilize the principe of an AND gate
    return accept_boundary * accept_sampler;
}


/* -------------------------------------------------------
   Sample 'nmoves' random moves
---------------------------------------------------------- */

void Sampler::sample(int nmoves)
{
    // declare variables
    int acceptance_counter = 0;

    // sample nmoves moves
    for (int i=0; i<nmoves; i++){
        // pick move type
        move_idx = rng->choice(system->moves_prob);
        ndrawn[move_idx] ++;
        Moves* move = system->moves[move_idx];

        //Moves* move = propose_move(box->moves, box->moves_prob);
        move->perform_move();
        if (!accept_move(move, system->temp, system->chempot)){
            move->reset();
        }
        else{
            naccepted[move_idx] ++;
            acceptance_counter ++;
        }
        for(Box* box : move->boxes){
            if (box->npar - 1 > box->nsystemsize.size()) {
                box->nsystemsize.resize(box->npar + 1);
            }
            box->nsystemsize[box->npar] ++;
        }
    }
    acceptance_ratio = (double)(acceptance_counter)/nmoves;
}
