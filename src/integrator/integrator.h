#pragma once
#include <iostream>
#include <armadillo>


using namespace std;
using namespace arma;


class Integrator
{
public:
    Integrator(class Box* box_in, double dt_in);
    virtual void next_step() = 0;
    virtual ~Integrator() = default;

protected:
    double dt;

    class Box* box = nullptr;
};
