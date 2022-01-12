#pragma once
#include <iostream>
#include <valarray>
#include <string>
#include <vector>


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
    std::vector<int> build_neigh_list(std::vector<class Particle>, int, double);
    void check_neighbors(int, std::vector<Particle>, unsigned int, std::vector<int>&,
                         std::vector<class Particle>, double);
    std::vector<int> detect_molecule(std::vector<class Particle>,
                                     std::vector<class Particle>,
                                     bool&, double);

    const double pi = 3.14159265358979323846;
    double du;

    class System* system = nullptr;
    class RandomNumberGenerator* rng = nullptr;
};
