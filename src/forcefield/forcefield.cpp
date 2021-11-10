#include "forcefield.h"
#include "../box.h"

ForceField::ForceField(class Box* box_in)
{
    box = box_in;
}

void ForceField::distance_matrix_element(const mat positions, const int i, const int j)
{
    /* Update element (i, j) of distance matrix and distance direction
     * cube. This corresponds to the distance between particles i and
     * j squared and the vector pointing from i to j, respectively.
     */
    rowvec dr = positions.row(j) - positions.row(i);
    distance_dir_cube.tube(i, j) = dr;
    distance_dir_cube.tube(j, i) = -dr;
    distance_mat(i, j) = sum(dr%dr);
    distance_mat(j, i) = distance_mat(i, j);
}


void ForceField::distance_matrix_cross(const mat positions, const int i)
{
    /* Update the distances between particle i and all other particles.
     * This corresponds to updating a cross in the distance matrix:
     * update row i and col i. Only upper triangular elements are filled.
     */
    for(int j=0; j<i; j++){
        distance_matrix_element(positions, i, j);
    }
    for(int j=i+1; j<box->npar; j++){
        distance_matrix_element(positions, j, i);
    }
}


void ForceField::add_distance_cross(const mat positions)
{
    /* Add row and column to distance matrix
     * when particle is added.
     */
    cout << "add_distance_cross1" << endl;
    distance_mat.insert_cols(box->npar-2, 1);
    cout << "add_distance_cross2" << endl;
    distance_mat.insert_rows(box->npar-2, 1);
    cout << "add_distance_cross3" << endl;
    distance_dir_cube.insert_cols(box->npar-2, 1);
    cout << "add_distance_cross4" << endl;
    distance_dir_cube.insert_rows(box->npar-2, 1);
    cout << "add_distance_cross5" << endl;
    for(int j=0; j<box->npar; j++){
        distance_matrix_element(positions, box->npar-1, j);
    }
}


void ForceField::rm_distance_cross(const int i){
    /* Remove row and column of distance matrix
     * when particle is removed.
     */
    distance_mat.shed_row(i);
    distance_mat.shed_col(i);
    distance_dir_cube.shed_row(i);
    distance_dir_cube.shed_col(i);
}


void ForceField::distance_matrix(const mat positions)
{
    /* Update all elements of distance matrix and distance
     * direction cube. Only upper triangular elements are filled.
     */
    distance_mat = zeros(box->npar, box->npar);
    distance_dir_cube = zeros(box->npar, box->npar, box->ndim);
    for(int i=0; i<box->npar; i++){
        for(int j=0; j<i; j++){
            distance_matrix_element(positions, i, j);
        }
    }
}


std::vector<int> ForceField::build_neigh_list(const int i, const double r_sqrd)
{
    /* Build neighbor list of particle i for a cutoff distance
     * (squared) r_sqrd.
     */
    std::vector<int> neigh_list;
    for(int j=0; j<i; j++){
        if(distance_mat(i, j) < r_sqrd){
            neigh_list.push_back(j);
        }
    }
    return neigh_list;
}

//ForceField::~ForceField() {}
