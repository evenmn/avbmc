#pragma once
#include <cmath>
#include <string>
#include <valarray>
#include "moves.h"


class TransMH : public Moves
{
public:
    TransMH(class System *, class Box *, double = 0.01, double = 0.01);
    void perform_move();
    double accept(double, double);
    void reset();
    void update_nsystemsize();
    std::string repr();

private:
    unsigned int i;
    double dx, Ddt;
    std::valarray<double> eps, df;
    class Box* box = nullptr;
};
