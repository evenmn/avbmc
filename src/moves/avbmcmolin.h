#pragma once
#include <cmath>
#include <string>
#include <vector>
#include <valarray>

#include "moves.h"


class AVBMCMolIn : virtual public Moves
{
public:
    AVBMCMolIn(class System *, class Box *, std::vector<class Particle>,
    double = 0.9, double = 1.5, double = 1.3, bool = true, bool = false);
    void perform_move() override;
    double accept(double, double) override;
    void reset() override;
    void update_size_histogram() override;
    std::string repr() override;

private:
    unsigned int detect_target_molecule(bool &);
    std::vector<Particle> create_molecule();
    std::valarray<double> insertion_position(unsigned int, bool = false); 

    bool detected_target, energy_bias, target_mol;
    unsigned int natom, neigh_id_above, neigh_id_below, neigh_id_inner;
    double r_below, r_above, r_inner, r_belowsq, r_abovesq, v_in, nmolavg, natom_inv;
    std::vector<unsigned int> npartype_old;
    std::vector<class Particle> particles, particles_old;
    class Box* box = nullptr;
};
