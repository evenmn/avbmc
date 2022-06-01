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

    bool detected_out;
    std::vector<class Particle> molecule_out;

private:
    unsigned int detect_target_molecule(bool &);
    std::vector<unsigned int> detect_deletion_molecule(unsigned int, bool &);


    unsigned int natom, neigh_id_above, neigh_id_inner, n_in;
    bool detected_target, energy_bias, target_mol;
    double v_in, nmolavg, r_inner, natom_inv;
    std::vector<class Particle> molecule;
    class Box* box = nullptr;
};
