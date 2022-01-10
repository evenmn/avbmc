#pragma once

#include <string>
#include <vector>
#include <valarray>


class ForceField
{
public:
    ForceField(class System *);

    // Declare pure virtual functions
    virtual void read_param_file(std::string) = 0;
    virtual void sort_params() = 0;
    double comp_energy_par(std::vector<class Particle>, int);
    virtual double comp_energy_par(std::vector<class Particle>, int, std::valarray<double> &, bool) = 0;
    virtual ~ForceField() = default;
    
    // Declare global functions
    std::vector<int> build_neigh_list(int, double);
    double norm(std::valarray<double>);

    // Store state properties to avoid unnecessary computations
    std::string label;
    std::string paramfile;
    std::vector<std::string> label1_vec, label2_vec;

protected:
    class System* system = nullptr;
};
