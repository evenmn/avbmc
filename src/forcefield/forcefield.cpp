#include <iostream>
#include <vector>
#include <valarray>

#include "forcefield.h"
#include "../system.h"
#include "../particle.h"

ForceField::ForceField(System* system_in)
{
    system = system_in;
}

/*
void ForceField::initialize()
{
    // Initialize matrices with correct dimensionality and compute
    // initial separation distance and force.
    //
    tmp_npar = box->npar;
    distance_mat = zeros(tmp_npar, tmp_npar);
    distance_dir_cube = zeros(tmp_npar, tmp_npar, box->ndim);
    poteng_mat = zeros(tmp_npar, tmp_npar);
    force_mat = zeros(tmp_npar, tmp_npar, box->ndim);
    tmp_to_state();
    //update_force_all(); 
}
*/
/*
void ForceField::update_distance_all()
{
    // Update all elements in distance matrix, potential energy matrix
    // and the force matrix. For Monte Carlo simulations,  this is 
    // necessary only when initializing the simulation. For molecular
    // dynamics simulations, this is done for every time step.
    //
    for(int i=0; i<box->npar; i++){
        for(int j=0; j<i; j++){
            std::valarray<double> dr = box->particles[j]->r - box->particles[i]->r;
            //distance_dir_cube.tube(i, j) = dr;
            //distance_dir_cube.tube(j, i) = -dr;
            //distance_mat(i, j) = sum(dr%dr);
            //distance_mat(j, i) = distance_mat(i, j);
        }
    }
}
*/
/*
void ForceField::update_distance_element(const mat positions, const int i, const int j)
{
    // Update element directly.
    //
    rowvec dr = positions.row(j) - positions.row(i);
    tmp_distance_dir_cube.tube(i, j) = dr;
    tmp_distance_dir_cube.tube(j, i) = -dr;
    tmp_distance_mat(i, j) = sum(dr%dr);
    tmp_distance_mat(j, i) = distance_mat(i, j);
}
*/
/*
void ForceField::update_distance_par(const mat positions, const int i)
{
    // Update the distances between particle i and all other particles.
    // This corresponds to updating a cross in the distance matrix:
    // update row i and col i.
    //

    cout << "update_distance_par1" << endl;
    for(int j=0; j<i; j++){
        cout << "update_distance_parx " << j << endl;
        rowvec dr = positions.row(j) - positions.row(i);
        cout << "update_distance_parxx " << dr << endl;
        tmp_distance_dir_cube.tube(i, j) = dr;
        cout << "update_distance_parxxx " << tmp_distance_dir_cube.n_rows << endl;
        tmp_distance_dir_cube.tube(j, i) = -dr;
        cout << "update_distance_parxxxx " << tmp_distance_mat.n_rows << endl;
        tmp_distance_mat(i, j) = sum(dr%dr);
        cout << "update_distance_parxxxxx " << tmp_distance_dir_cube.n_slices << endl;
        tmp_distance_mat(j, i) = tmp_distance_mat(i, j);
    }
    cout << "update_distance_par2" << endl;
    for(int j=i+1; j<tmp_npar; j++){
        rowvec dr = positions.row(j) - positions.row(i);
        tmp_distance_dir_cube.tube(i, j) = dr;
        tmp_distance_dir_cube.tube(j, i) = -dr;
        tmp_distance_mat(i, j) = sum(dr%dr);
        tmp_distance_mat(j, i) = tmp_distance_mat(i, j);
    }
    cout << "update_distance_par3" << endl;
}
*/
/*
void ForceField::add_distance_par(const mat positions, const std::string chem_symbol, const int type, const double mass)
{
    // Add row and column to distance matrix
    // when particle is added.
    //
    tmp_npar ++;
    tmp_to_state();

    tmp_distance_mat.insert_cols(box->npar-1, 1);
    tmp_distance_mat.insert_rows(box->npar-1, 1);
    tmp_distance_dir_cube.insert_cols(box->npar-1, 1);
    tmp_distance_dir_cube.insert_rows(box->npar-1, 1);
    tmp_poteng_mat.insert_cols(box->npar-1, 1);
    tmp_poteng_mat.insert_rows(box->npar-1, 1);
    tmp_force_mat.insert_cols(box->npar-1, 1);
    tmp_force_mat.insert_rows(box->npar-1, 1);
    //tmp_chem_symbols.push_back(chem_symbol);
    //tmp_particle_types.push_back(type);
    //tmp_particle_masses.push_back(mass);
    update_distance_par(positions, box->npar);
}
*/
/*
void ForceField::rm_distance_par(const int i){
    // Remove row and column of temporary distance matrix
    // when particle is removed, but before it is accepted.
    //
    tmp_npar --;
    tmp_to_state();

    tmp_distance_mat.shed_row(i);
    tmp_distance_mat.shed_col(i);
    tmp_distance_dir_cube.shed_row(i);
    tmp_distance_dir_cube.shed_col(i);
    tmp_poteng_mat.shed_row(i);
    tmp_poteng_mat.shed_col(i);
    tmp_force_mat.shed_row(i);
    tmp_force_mat.shed_col(i);
    //tmp_chem_symbols.erase(tmp_chem_symbols.begin() + i);
    //tmp_particle_masses.erase(tmp_particle_masses.begin() + i);
    //tmp_particle_types.erase(tmp_particle_types.begin() + i);
    //update_distance_par(positions, box->npar);
}
*/
/*
void ForceField::tmp_to_state()
{
    // Set temporary matrices equal to the state matrices
    //
    tmp_distance_mat = distance_mat;
    tmp_distance_dir_cube = distance_dir_cube;
    tmp_poteng_mat = poteng_mat;
    tmp_force_mat = force_mat;
    //tmp_chem_symbols = box->chem_symbols;
    //tmp_particle_types = box->particle_types;
    //tmp_particle_masses = box->particle_masses;
}
*/
/*
void ForceField::state_to_tmp()
{
    // Set state matrices equal to the temporary matrices
    //
    distance_mat = tmp_distance_mat;
    distance_dir_cube = tmp_distance_dir_cube;
    poteng_mat = tmp_poteng_mat;
    force_mat = tmp_force_mat;
    //box->chem_symbols = tmp_chem_symbols;
    //box->particle_types = tmp_particle_types;
    //box->particle_masses = tmp_particle_masses;
}
*/
/*
void ForceField::distance_matrix(const mat positions)
{
    // Update all elements of distance matrix and distance
    // direction cube. Only upper triangular elements are filled.
    //
    distance_mat = zeros(box->npar, box->npar);
    distance_dir_cube = zeros(box->npar, box->npar, box->ndim);
    for(int i=0; i<box->npar; i++){
        for(int j=0; j<i; j++){
            distance_matrix_element(positions, i, j);
        }
    }
}
*/

/* -------------------------------------------------------------
   Build neighbor list of particle 'i' with maximum neighbor
   distance squared 'rsq'
---------------------------------------------------------------- */
/*
std::vector<int> ForceField::build_neigh_list(const int i, const double rsq)
{
    double rijsq;
    std::valarray<double> ri = box->particles[i]->r;
    std::vector<int> neigh_list;
    for(int j=0; j<i; j++){
        rijsq = std::pow(box->particles[j]->r - posi, 2).sum();
        if(sdist < r_sqrd){
            neigh_list.push_back(j);
        }
    }
    for(int j=i+1; j<box->npar; j++){
        sdist = std::pow(box->particles[j]->r - posi, 2).sum();
        if(sdist < r_sqrd){
            neigh_list.push_back(j);
        }
    }
    return neigh_list;
}
*/

/* -------------------------------------------------------------
   Compute the squared norm of a valarray 'array'
---------------------------------------------------------------- */

double ForceField::norm(std::valarray<double> array)
{
    double normsq = 0.;
    for (unsigned int i=0; i < array.size(); i++)
    {
        normsq += array[i] * array[i];
    }
    return normsq;
}


//ForceField::~ForceField() {}
