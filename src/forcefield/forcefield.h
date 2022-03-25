#pragma once

#include <map>
#include <string>
#include <vector>
#include <valarray>


class ForceField
{
public:
    ForceField(class Box *);

    // Declare pure virtual functions
    virtual void read_param_file(std::string) = 0;
    virtual void allocate_memory() = 0;
    virtual void free_memory() = 0;
    virtual double comp_energy_par_neigh0_eng0(int, std::valarray<double> &, bool) = 0;
    virtual double comp_energy_par_neigh1_eng0(int, std::valarray<double> &, bool) = 0;
    virtual double comp_energy_par_neigh1_eng1(int, std::valarray<double> &, bool) = 0;
    virtual ~ForceField() = default;
    
    // Declare global functions
    void initialize(), set(), reset();
    double comp_energy_par_force0(int);
    double comp_energy_par_force1(int, std::valarray<double> &);
    double (ForceField::*comp_energy_par)(int, std::valarray<double> &, bool) = nullptr;

    // Store state properties to avoid unnecessary computations
    std::map<std::string, unsigned int> label2type;
    unsigned int ntype;
    std::string label, paramfile;
    std::vector<std::string> label1_vec, label2_vec, unique_labels;
    double temp_scale = 1.;
    std::vector<double> poteng_vec, poteng_vec_old;
    std::vector<std::valarray<double> > force_vec, force_vec_old;
    std::vector<std::vector<double> > poteng_mat, poteng_mat_old;
    std::vector<std::vector<std::valarray<double > > > force_cube, force_cube_old;

protected:
    double norm(std::valarray<double>);
    void create_label_mapping();
    class Box* box = nullptr;

    unsigned int nline;
};
