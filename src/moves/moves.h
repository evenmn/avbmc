#pragma once
#include <iostream>
#include <valarray>
#include <vector>


class Moves
{
public:
    Moves(class Box *);

    // declare pure virtual functions
    virtual void perform_move() = 0;
    virtual double accept(double, double) = 0;
    virtual void reset() = 0;
    virtual ~Moves() = default;

protected:
    double pi = 3.14159265358979323846;
    double du;

    class Box* box = nullptr;
    class RandomNumberGenerator* rng = nullptr;
};
