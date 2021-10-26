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

mat Sampler::perform_move(class Moves* move, const int j)
{
    /* Perform move 'move' on particle j
     */
    return move->perform_move(j);
}

bool Sampler::accept_move(class Moves* move, double temp, double chempot)
{
    /* Decide if move should be accepted
     * or rejected
     */
    double p = get_acceptance_prob(move, temp, chempot);
    //cout << p << endl;
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
        int j = rng->next_int(box->npar); // might move back to moves
        move->perform_move(j);
        accepted = accept_move(move, temp, chempot);
        if(accepted){
            move->update_box(j);
            acceptance_counter ++;
        }
    }
    acceptance_ratio = (double)(acceptance_counter)/nmoves;
}
