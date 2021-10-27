#pragma once
#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;


class Boundary
{
public:
    Boundary(class Box* box_in);
    virtual bool correct_position(mat &pos) = 0;
    virtual bool correct_velocity(mat &vel) = 0;   // needed for reflective boundaries
    virtual bool correct_distance(mat &dist) = 0;  // needed for periodic boundaries
    virtual double comp_volume() = 0;
    virtual ~Boundary() = default;

protected:
    class Box* box = nullptr;
};
