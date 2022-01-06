#pragma once

#include <string>
#include <vector>
#include <valarray>
#include <memory>


class ForceField
{
public:
    ForceField(class System *);

    // Declare pure virtual functions
    virtual void read_param_file(std::string) = 0;
    virtual void sort_params() = 0;
    //virtual void build_neigh_lists() = 0;
    //virtual double eval_acc_element(const mat positions, const int i, const int j, rowvec& acc, const bool comp_energy=false) = 0;
    //virtual double eval_acc_par(const mat positions, const int i, rowvec& acc, const bool comp_energy=false) = 0;
    //virtual double eval_acc(const mat positions, mat& accs, vec& potengs, const bool comp_energy=false) = 0;
    //virtual double comp_force_par(const rowvec pos, rowvec& acc) = 0;
    //virtual double update_force_par(const mat positions, const int i) = 0;
    //virtual double update_force_all() = 0;
    //virtual double comp_energy_mol(std::vector<class Particle *>, class Molecule *) = 0;
    //virtual double comp_energy_par(std::vector<class Particle *>, int) = 0;
    virtual double comp_energy_mol(std::vector<class Particle>, class Molecule *) = 0;
    virtual double comp_energy_par(std::vector<class Particle>, int) = 0;
    virtual ~ForceField() = default;
    
    // Declare global functions
    //void distance_matrix_element(const mat positions, const int i, const int j);
    //void distance_matrix_cross(const mat positions, const int i);
    //void distance_matrix(const mat positions);
    //void initialize();
    //void update_distance_par(const mat positions, const int i);
    //void update_distance_all();
    //void add_distance_par(const mat positions, const std::string label, const int type, const double mass);
    //void rm_distance_par(const int i);
    //void tmp_to_state();
    //void state_to_tmp();
    std::vector<int> build_neigh_list(int, double);
    double norm(std::valarray<double>);

    // Store state properties to avoid unnecessary computations
    // Matrices have dimensionality (npar, npar)
    //mat distance_mat, poteng_mat;
    //cube distance_dir_cube, force_mat;

    // Store temporary state properties. This is convenient for
    // Monte Carlo simulations when we do not know if the 
    // state will be accepted
    //mat tmp_distance_mat, tmp_poteng_mat, tmp_accelerations;
    //cube tmp_distance_dir_cube, tmp_force_mat;

    //std::vector<std::string> tmp_chem_symbols;
    //std::vector<int> tmp_particle_types;
    //std::vector<double> tmp_particle_masses;
    
    std::string label;
    std::string paramfile;
    std::vector<std::string> label1_vec, label2_vec;

protected:
    class System* system = nullptr;
};
