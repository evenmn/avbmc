#include <iostream>
#include <string>
#include <armadillo>

using namespace arma;

class Vashishta : public ForceField
{
public:
    Vashishta(string params);
    double eval_acc_par(const vec positions, const int i, vec& acc, const bool comp_energy=false) = 0;
    double eval_acc(const vec positions, vec& acc, const bool comp_energy=false) = 0;
    
private:
    string label = "Vashishta";
};
