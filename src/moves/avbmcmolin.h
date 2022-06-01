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
    void insert(std::vector<class Particle>);

private:
    unsigned int detect_target_molecule(bool &);
    std::vector<Particle> create_molecule();

    bool detected_target, energy_bias, target_mol;
    unsigned int natom, neigh_id_above, neigh_id_below, neigh_id_inner;
    double r_inner, v_in, nmolavg, natom_inv;
    std::vector<class Particle> particles;
    class Box* box = nullptr;
};
