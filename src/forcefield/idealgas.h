#pragma once

#include <string>
#include <vector>
#include <valarray>

#include "forcefield.h"


class IdealGas : public ForceField
{
public:
    IdealGas(class Box *, std::vector<std::string>);
    
    double comp_energy_par_neigh0_eng0(int, std::valarray<double> &, bool);
    double comp_energy_par_neigh1_eng0(int, std::valarray<double> &, bool);
    double comp_energy_par_neigh1_eng1(int, std::valarray<double> &, bool);
};
