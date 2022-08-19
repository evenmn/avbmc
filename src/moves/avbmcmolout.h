/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.
---------------------------------------------------------------------------- */

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
    unsigned int nrejecttarget, nrejectout;
    std::vector<class Particle> molecule_out;

private:
    void check_neigh_recu(int, unsigned int, std::vector<unsigned int> &, std::vector<std::vector<unsigned int> >);
    unsigned int detect_target_molecule(bool &);
    std::vector<unsigned int> detect_deletion_molecule(std::vector<unsigned int>, bool &);


    //std::vector<std::vector<unsigned int> > *neigh_list_above, *neigh_list_inner;
    unsigned int natom, neigh_id_above, neigh_id_inner, n_in, center_type;
    bool detected_target, energy_bias, target_mol;
    double v_in, nmolavg, r_inner, natom_inv;
    std::vector<class Particle> molecule;
    class Box* box = nullptr;
};
