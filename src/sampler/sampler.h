#pragma once
#include <iostream>
#include <functional>
#include <vector>

using namespace std;


class Sampler
{
public:
    Sampler(class Box* box_in);

    // declare pure virtual functions
    virtual double w(int npar) = 0;
    virtual ~Sampler() = default;
    
    // declare public functions
    class Moves* propose_move(std::vector<class Moves*> moves, std::vector<double> move_probs);
    bool accept_move(class Moves* move, const double temp, const double chempot);
    void sample(int nmoves);

    // declare public variables
    int move_idx;
    double acceptance_ratio;

protected:
    class Box* box = nullptr;
    class RandomNumberGenerator* rng = nullptr;
};
