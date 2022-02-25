#pragma once

#include <map>
#include <string>
#include <vector>
#include <valarray>


class ForceField
{
public:
    ForceField(class System *);

    // Declare pure virtual functions
    virtual void read_param_file(std::string) = 0;
    virtual void allocate_memory() = 0;
    virtual void free_memory() = 0;
    virtual void sort_params() = 0;
    virtual double comp_energy_par(std::vector<class Particle>, int, std::valarray<double> &, bool) = 0;
    virtual ~ForceField() = default;
    
    // Declare global functions
    double comp_energy_par(std::vector<class Particle>, int);
    std::vector<int> build_neigh_list(int, double);

    // Store state properties to avoid unnecessary computations
    std::map<std::string, unsigned int> label2type;
    unsigned int ntype;
    std::string label, paramfile;
    std::vector<std::string> label1_vec, label2_vec, unique_labels;
    double temp_scale = 1.;

protected:
    double norm(std::valarray<double>);
    void create_label_mapping();
    class System* system = nullptr;

    unsigned int nline;
};
