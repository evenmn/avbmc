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
   Decide if move should be accepted or rejected. This
   is done in two steps:
     1. Check if boundary condition is satisfied
     2. Check if acceptance criterion is satisfied
--------------------------------------------------------- */

bool Sampler::accept_move(Moves* move, double temp, double chempot)
{
    std::cout << "accept_move1" << std::endl;
    // check boundary condition
    bool accept_boundary = true;
    std::cout << "accept_move2" << std::endl;
    for(Box* box : move->boxes){
        accept_boundary *= box->boundary->correct_position();
    }
    std::cout << "accept_move3" << std::endl;

    // check move acceptance probability
    bool accept_sampler = move->accept(temp, chempot) > rng->next_double();
    std::cout << "accept_move4" << std::endl;

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
        std::cout << "sample1" << std::endl;
        move_idx = rng->choice(system->moves_prob);
        std::cout << "sample2" << std::endl;
        Moves* move = system->moves[move_idx];
        std::cout << move->label << std::endl;
        std::cout << "sample3" << std::endl;
        move->ndrawn ++;
        std::cout << "sample4" << std::endl;

        move->perform_move();
        std::cout << "sample5" << std::endl;
        if (!accept_move(move, system->temp, system->chempot)){
            move->reset();
        }
        else{
            move->naccept ++;
            acceptance_counter ++;
        }
        std::cout << "sample6" << std::endl;
        for(Box* box : move->boxes){
            if (box->npar - 1 > box->nsystemsize.size()) {
                box->nsystemsize.resize(box->npar + 1);
            }
            box->nsystemsize[box->npar] ++;
        }
        std::cout << "sample7" << std::endl;
    }
    acceptance_ratio = (double)(acceptance_counter)/nmoves;
}
