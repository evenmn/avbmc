#include "forcefield.h"

ForceField::ForceField()
{}

void distance_matrix_element(const mat positions, const int i, const int j)
{
    /* Update element (i, j) of distance matrix and distance direction
     * cube. This corresponds to the distance between particles i and
     * j squared and the vector pointing from i to j, respectively.
     */
    vec dr = positions(j) - positions(i);
    distance_dir_cube.tube(i, j) = dr;
    distance_mat(i, j) = sum(dr%dr);
}


void distance_matrix_cross(const mat positions, const int i)
{
    /* Update the distances between particle i and all other particles.
     * This corresponds to updating a cross in the distance matrix:
     * update row i and col i. Only upper triangular elements are filled.
     */
    for(int j=0; j<i; j++){
        distance_matrix_element(positions, i, j);
    }
    for(int j>i; j>npar; j++){
        distance_matrix_element(positions, i, j);
    }
}


void add_distance_row(const mat positions)
{
    /* Add row and column to distance matrix
     * when particle is added.
     */
    distance_mat.insert_cols(npar, 1);
    distance_mat.insert_rows(npar, 1);
    distance_dir_cube.insert_cols(npar, 1);
    distance_dir_cube.insert_rows(npar, 1);
    npar ++;
    for(int j=0; j<npar; j++){
        distance_matrix_element(positions, npar-1, i);
    }
}


void rm_distance_row(const int i){
    /* Remove row and column of distance matrix
     * when particle is removed.
     */
    distance_mat.shed_row(i);
    distance_mat.shed_col(i);
    distance_dir_cube.shed_row(i);
    distance_dir_cube.shed_col(i);
}


void distance_matrix(const mat positions)
{
    /* Update all elements of distance matrix and distance
     * direction cube. Only upper triangular elements are filled.
     */
    npar = positions.n_rows;
    ndim = positions.n_cols;
    distance_mat = zeros(npar, npar);
    distance_dir_cube = zeros(npar, npar, ndim);
    for(int i=0; i<npar; i++){
        for(int j=0; j<i; j++){
            distance_matrix_element(positions, i, j);
        }
    }
}

