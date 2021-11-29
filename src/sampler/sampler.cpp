#include "sampler.h"
#include "../box.h"
#include "../moves/moves.h"

Sampler::Sampler(Box* box_in)
{
    box = box_in;
    rng = box->rng;

    acceptance_ratio = 0.;
}

class Moves* Sampler::propose_move(std::vector<Moves*> moves, std::vector<double> moves_prob)
{
    /* Given a vector of moves and a vector of corresponding
     * probabilities, returning move
     */
    move_idx = rng->choice(moves_prob);
    return moves[move_idx];
}

bool Sampler::accept_move(Moves* move, double temp, double chempot)
{
    /* Decide if move should be accepted or rejected.
     * This is done in two steps:
     *  1. Check if boundary condition is satisfied
     *  2. Check if sampler condition is satisfied
     */

    // check boundary condition
    bool accept_boundary = box->boundary->correct_position();

    // check move acceptance probability
    bool accept_sampler = move->accept(temp, chempot) > rng->next_double();

    // utilize the principe of an AND gate
    return accept_boundary * accept_sampler;
}

void Sampler::sample(int nmoves)
{
    /* Sample 
     */

    // declare variables
    int acceptance_counter = 0;

    // sample nmoves moves
    for(int i=0; i<nmoves; i++){
        Moves* move = propose_move(box->moves, box->moves_prob);
        move->perform_move();
        if(!accept_move(move, box->temp, box->chempot)){
            move->reset();
        }
        else{
            acceptance_counter ++;
        }
    }
    acceptance_ratio = (double)(acceptance_counter)/nmoves;
}
