#pragma once
#include <vector>


class Sampler
{
public:
    Sampler(class Box*);

    // declare pure virtual functions
    virtual double w(int) = 0;
    virtual ~Sampler() = default;
    
    // declare public functions
    void sample(int);

    // declare public variables
    int move_idx;
    double acceptance_ratio;

private:
    // declare private functions
    class Moves* propose_move(std::vector<class Moves*>, std::vector<double>);
    bool accept_move(class Moves*, double, double);

protected:
    class Box* box = nullptr;
    class RandomNumberGenerator* rng = nullptr;
};
