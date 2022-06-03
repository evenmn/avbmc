/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.
---------------------------------------------------------------------------- */

#pragma once
#include <string>
#include <vector>
#include <valarray>

#include "forcefield.h"


class IdealGas : public ForceField
{
public:
    IdealGas(class Box *, const std::vector<std::string> &);
    
    double comp_energy_par_neigh0_eng0(unsigned int,
        std::valarray<double> &, bool) override;
    double comp_energy_par_neigh1_eng0(unsigned int,
        std::valarray<double> &, bool) override;
    double comp_energy_par_neigh1_eng1(unsigned int,
        std::valarray<double> &, bool) override;
};
