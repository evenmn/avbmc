#include "sampler.h"
#include "../box.h"
#include "../moves/moves.h"

Sampler::Sampler(class Box* box_in)
{
    box = box_in;
    rng = box->rng;

    acceptance_ratio = 0.;
    move_idx = NAN;
}

class Moves* Sampler::propose_move(vector<class Moves*> moves, std::vector<double> moves_prob)
{
    /* Given a vector of moves and a vector of corresponding
     * probabilities, returning move
     */
    move_idx = rng->choice(moves_prob);
    return moves[move_idx];
}

bool Sampler::accept_move(class Moves* move, double temp, double chempot)
{
    /* Decide if move should be accepted or rejected.
     * This is done in two steps:
     *  1. Check if boundary condition is satisfied
     *  2. Check if sampler condition is satisfied
     */

    // check boundary condition
    bool accept_boundary = box->boundary->correct_position();

    // check sampler condition
    double p = get_acceptance_prob(move, temp, chempot);
    std::cout << p << std::endl;
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
        move->perform_move();
        accepted = accept_move(move, box->temp, box->chempot);
        if(!accepted){
            move->reset();
        }
        else{
            acceptance_counter ++;
        }
    }
    acceptance_ratio = (double)(acceptance_counter)/nmoves;
}
