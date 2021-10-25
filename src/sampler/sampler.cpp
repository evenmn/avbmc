#include "sampler.h"
#include "../box.h"

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

vec Sampler::perform_move(class Moves* move)
{
    /* Propose move
     */
    return move->perform_move();
}

bool Sampler::accept_move(class Moves* move, double temp, double chempot)
{
    /* Decide if move should be accepted
     * or rejected
     */
    double p = get_acceptance_prob(move, temp, chempot);
    if(p > rng->next_double()){
        return true;
    }
    return false;
}

void Sampler::sample(int nmoves)
{
    /* Sample 
     */

    // declare variables
    bool accepted;
    int acceptance_counter = 0;

    double temp = box->temp;
    double chempot = box->chempot;

    class Moves* move = nullptr;
    vector<class Moves*> moves = box->moves;
    vector<double> moves_prob = box->moves_prob;

    // sample nmoves moves
    for(int i=0; i<nmoves; i++){
        move = propose_move(moves, moves_prob);
        vec pos_new = perform_move(move);
        accepted = accept_move(move, temp, chempot);
        if(accepted){
            box->positions = pos_new;
            box->npar = pos_new.n_rows;
            box->poteng += du;
            acceptance_counter ++;
        }
    }
    acceptance_ratio = (double)(acceptance_counter/nmoves);
}
