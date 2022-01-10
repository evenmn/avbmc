#pragma once
#include <iostream>
#include <valarray>
#include <string>
#include <vector>
#include <memory>


class Moves
{
public:
    Moves(class System *);

    // declare pure virtual functions
    virtual void perform_move() = 0;
    virtual double accept(double, double) = 0;
    virtual void reset() = 0;
    virtual void update_nsystemsize() = 0;
    virtual std::string repr() = 0;
    virtual ~Moves() = default;

    int ndrawn, naccept;
    std::string label;

protected:
    std::vector<std::valarray<double> > rotate_molecule(std::vector<std::valarray<double> >);
    double norm(std::valarray<double>);

    const double pi = 3.14159265358979323846;
    double du;

    class System* system = nullptr;
    class RandomNumberGenerator* rng = nullptr;
};
