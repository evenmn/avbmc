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
    distance_mat(i, j) = sum(dr%dr);
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
    for(int j=i+1; j>box->npar; j++){
        distance_matrix_element(positions, i, j);
    }
}


void ForceField::add_distance_cross(const mat positions)
{
    /* Add row and column to distance matrix
     * when particle is added.
     */
    distance_mat.insert_cols(box->npar, 1);
    distance_mat.insert_rows(box->npar, 1);
    distance_dir_cube.insert_cols(box->npar, 1);
    distance_dir_cube.insert_rows(box->npar, 1);
    box->npar ++;
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

//ForceField::~ForceField() {}
