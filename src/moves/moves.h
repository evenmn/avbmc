#pragma once
#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

class Moves
{
public:
    Moves(class Box* box_in);

    // declare pure virtual functions
    virtual mat perform_move(const int i) = 0;
    virtual double accept() = 0;
    virtual void update_box(const int i) = 0;
    virtual ~Moves() = default;

    // declare public functions
    bool check_stillinger();
    double get_cluster_volume();


protected:
    class Box* box = nullptr;
    class RandomNumberGenerator* rng = nullptr;
};
