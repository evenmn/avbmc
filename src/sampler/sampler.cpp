#include "sampler.h"
#include "../box.h"
#include "../moves/moves.h"

Sampler::Sampler(class Box* box_in)
{
    box = box_in;
    rng = box->rng;
}

class Moves* Sampler::propose_move(vector<class Moves*> moves, vector<double> moves_prob)
{
    /* Given a vector of moves and a vector of corresponding
     * probabilities, returning move
     */
    int i = rng->choice(moves_prob);
    return moves[i];
}
/*
mat Sampler::perform_move(class Moves* move, const int j)
{
    // Perform move 'move' on particle j
    return move->perform_move(j);
}
*/

bool Sampler::accept_move(class Moves* move, double temp, double chempot)
{
    /* Decide if move should be accepted or rejected.
     * This is done in two steps:
     *  1. Check if boundary condition is satisfied
     *  2. Check if sampler condition is satisfied
     */

    // check boundary condition
    bool accept_boundary = box->boundary->correct_position(box->positions);

    // check sampler condition
    double p = get_acceptance_prob(move, temp, chempot);
    bool accept_sampler = p > rng->next_double();

    // utilize the principe of an AND gate
    return accept_boundary * accept_sampler;
}

void Sampler::sample(int nmoves)
{
    /* Sample 
     */

    // declare variables
    bool accepted;
    int acceptance_counter = 0;
    class Moves* move = nullptr;

    // sample nmoves moves
    for(int i=0; i<nmoves; i++){
        move = propose_move(box->moves, box->moves_prob);
        int j = rng->next_int(box->npar); // might move back to moves
        move->perform_move(j);
        accepted = accept_move(move, box->temp, box->chempot);
        if(accepted){
            move->update_box(j);
            acceptance_counter ++;
        }
    }
    acceptance_ratio = (double)(acceptance_counter)/nmoves;
}
