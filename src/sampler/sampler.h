#pragma once
#include <iostream>
#include <armadillo>

//#include "../moves/moves.h"

using namespace std;
using namespace arma;


class Sampler
{
public:
    Sampler(class Box* box_in);

    // declare pure virtual functions
    virtual double get_acceptance_prob(class Moves* move, double temp, double chempot) = 0;
    virtual ~Sampler() = default;
    
    // declare public functions
    class Moves* propose_move(vector<class Moves*> moves, vector<double> move_probs);
    mat perform_move(class Moves* move, const int j);
    bool accept_move(class Moves* move, const double temp, const double chempot);
    void sample(int nmoves);

    // declare public variables
    rowvec da;
    double du, acceptance_ratio;

protected:
    class Box* box = nullptr;
    class RandomNumberGenerator* rng = nullptr;
};
