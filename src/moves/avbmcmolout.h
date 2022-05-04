#pragma once
#include <string>
#include <vector>

#include "moves.h"


class AVBMCMolOut : virtual public Moves
{
public:
    AVBMCMolOut(class System *, class Box *, std::vector<Particle>,
                   double = 3.0, double = 2.0, bool = true, bool = false);
    void perform_move() override;
    double accept(double, double) override;
    void reset() override;
    void update_size_histogram() override;
    std::string repr() override;

private:
    unsigned int natom, neigh_id_above, neigh_id_inner;
    bool reject_move, energy_bias, target_mol;
    double r_above, r_abovesq, v_in, nmolavg, r_inner, natom_inv;
    std::vector<unsigned int> npartype_old;
    std::vector<class Particle> particles_old, molecule;
    class Box* box = nullptr;
};
