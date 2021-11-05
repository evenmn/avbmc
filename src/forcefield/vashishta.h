#pragma once
#include "forcefield.h"

class Vashishta : public ForceField
{
public:
    Vashishta(Box box, string params);
    void read_param_file(const string params);
    double eval_acc_par(const vec positions, const int i, vec& acc, const bool comp_energy=false) = 0;
    double eval_acc(const vec positions, vec& acc, const bool comp_energy=false) = 0;
    
private:
    string label = "Vashishta";
};
