#include <iostream>
#include <string>
#include <vector>
#include <armadillo>

using namespace arma;

class LennardJones : public ForceField
{
public:
    LennardJones(string params);
    void read_param_file(string params);
    double eval_force_par(const mat positions, const int i, vec& acc, const bool comp_energy=false) = 0;
    double eval_force(const mat positions, vec& acc, const bool comp_energy=false) = 0;

private:
    double sigma, epsilon, rc, rc_sqrd;
    field<vector<int>> neigh_lists();
    string label = "Lennard Jones";
};
