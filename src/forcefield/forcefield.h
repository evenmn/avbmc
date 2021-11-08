#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <armadillo>

using namespace std;
using namespace arma;


class ForceField
{
public:
    ForceField(class Box* box_in);

    // Declare pure virtual functions
    virtual void read_param_file(const string params) = 0;
    virtual void sort_params() = 0;
    virtual void build_neigh_lists(const mat positions) = 0;
    virtual double eval_acc_element(const mat positions, const int i, const int j, rowvec& acc, const bool comp_energy=false) = 0;
    virtual double eval_acc_par(const mat positions, const int i, rowvec& acc, const bool comp_energy=false) = 0;
    virtual double eval_acc(const mat positions, mat& accs, vec& potengs, const bool comp_energy=false) = 0;
    virtual ~ForceField() = default;
    
    // Declare global functions
    void distance_matrix_element(const mat positions, const int i, const int j);
    void distance_matrix_cross(const mat positions, const int i);
    void distance_matrix(const mat positions);
    void add_distance_cross(const mat positions);
    void rm_distance_cross(const int i);

    mat distance_mat;
    cube distance_dir_cube;

protected:
    class Box* box = nullptr;
};
