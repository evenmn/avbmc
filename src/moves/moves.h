/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.
---------------------------------------------------------------------------- */

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
    virtual void update_size_histogram() = 0;
    virtual std::string repr() = 0;
    virtual ~Moves() = default;

    int ndrawn, naccept;
    double cum_time;
    std::string label;

protected:
    std::vector<class Particle> rotate_molecule(std::vector<class Particle>);
    std::valarray<double> insertion_position(bool = false);
    double norm(std::valarray<double>);
    std::vector<unsigned int> build_neigh_list(std::vector<class Particle>,
        unsigned int, double);
    void check_neigh_recu(int, std::vector<class Particle>, unsigned int,
        std::vector<unsigned int> &, std::vector<class Particle>, double);
    void check_neigh_recu(int, std::vector<class Particle>, unsigned int,
        std::vector<unsigned int> &, std::vector<std::vector<unsigned int> >);
    void check_neigh_recu(int, std::vector<class Particle>,
        std::vector<class Particle>, unsigned int, std::vector<unsigned int> &,
        std::vector<std::vector<unsigned int> >);
    std::vector<unsigned int> detect_molecule(std::vector<class Particle>,
        std::vector<class Particle>, bool &, double);
    std::vector<unsigned int> detect_molecule(std::vector<std::vector<unsigned int> >,
        std::vector<class Particle>, bool &);
    std::vector<unsigned int> detect_molecule(std::vector<std::vector<unsigned int> >,
        const std::vector<class Particle> &, std::vector<class Particle>, bool &);

    const double pi = 3.14159265358979323846;
    double du, r_above, r_abovesq, r_below, r_belowsq;

    class System* system = nullptr;
    class RandomNumberGenerator* rng = nullptr;
};
