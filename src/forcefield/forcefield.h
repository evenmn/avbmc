#include <iostream>
#include <string>
#include <armadillo>

using namespace arma;

class ForceField
{
public:
    ForceField();
    virtual double eval_force_par(const mat positions, const int i, vec& acc, const bool comp_energy=false) = 0;
    virtual double eval_force(const mat positions, vec& acc, const bool comp_energy=false) = 0;
        
    void distance_matrix_element(const mat positions, const int i, const int j);
    void distance_matrix_cross(const mat positions, const int i);
    void distance_matrix(const mat positions);
    void add_distance_col(const mat positions);
    void rm_distance_col(const int i);

    mat distance_mat;
    cube distance_dir_cube;

protected:
    int npar;
    int ndim;
};
