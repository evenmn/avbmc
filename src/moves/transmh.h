#pragma once
#include <cmath>
#include <string>
#include <valarray>
#include "moves.h"


class TransMH : public Moves
{
public:
    TransMH(class System *, class Box *, double = 0.01, double = 0.01);
    void perform_move() override;
    double accept(double, double) override;
    void reset() override;
    void update_size_histogram() override;
    std::string repr() override;

private:
    unsigned int i;
    double dx, Ddt;
    std::valarray<double> eps, df;
    class Box* box = nullptr;
};
