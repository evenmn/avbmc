#pragma once
#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

class Moves
{
public:
    Moves(class Box* box_in);
    virtual mat perform_move() = 0;
    virtual double accept() = 0;
    virtual ~Moves() = default;

    bool check_stillinger();
    double get_cluster_volume();

protected:
    class Box* box = nullptr;
    class RandomNumberGenerator* rng = nullptr;
};
