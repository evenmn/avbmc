#include "lennardjones.h"
#include "../box.h"

LennardJones::LennardJones(class Box* box_in, string params)
    : ForceField(box_in)
{
    read_param_file(params);
}

void LennardJones::read_param_file(string params)
{
    /* Read parameter file and store parameters
     * globally. It takes the following form:
     * <type1> <type2> <sigma> <epsilon> <rc>
     * Parameters are stored in matrices of 
     * size (num_types x num_types)
     */

    // replace this with a parser
    vector<int> type1_vec = {0};
    vector<int> type2_vec = {0};
    vector<double> sigma_vec = {1.};
    vector<double> epsilon_vec = {1.};
    vector<double> rc_vec = {3.};

    int max_type1 = *max_element(type1_vec.begin(), type1_vec.end());
    int max_type2 = *max_element(type2_vec.begin(), type2_vec.end());
    int max_type = max(max_type1, max_type2) + 1;

    sigma_mat.zeros(max_type, max_type);
    epsilon_mat.zeros(max_type, max_type);
    rc_sqrd_mat.zeros(max_type, max_type);

    double rc;
    int type1, type2;
    for(int i=0; i<type1_vec.size(); i++){
        type1 = type1_vec[i];
        type2 = type2_vec[i];
        sigma_mat(type1, type2) = sigma_vec[i];
        epsilon_mat(type1, type2) = epsilon_vec[i];
        rc = rc_vec[i]; 
        rc_sqrd_mat(type1, type2) = rc * rc;
    }
}


vector<int> LennardJones::build_neigh_list(const mat positions, const int i)
{
    /* Build neighbor list for a particle i
     */
    vector<int> neigh_list;
    int typei = box->types(i);
    for(int j=0; j<i; j++){
        int typej = box->types(j);
        if(distance_mat(i, j) < rc_sqrd_mat(typei, typej)){
            neigh_list.push_back(j);
        }
    }
    return neigh_list;
}


void LennardJones::build_neigh_lists(const mat positions)
{
    /* Build neigh lists for all particles
     */
    for(int i=0; i<box->npar; i++){
        neigh_lists.push_back(build_neigh_list(positions, i));
    }
}


double LennardJones::eval_acc_par(const mat positions, const int i, rowvec& acc, const bool comp_energy)
{
    /* Evaluate force acting on particle i. This is needed only when
     * particle i is moved, and therefore we also have to update the
     * distance between particle i and all other particles.
     */
    // Update distance
    //ForceField::distance_matrix_cross(positions, i);
    //build_neigh_list(positions, i);
    
    // Get type of particle i
    int typei, typej;
    typei = box->types(i);

    acc.zeros(box->ndim);
    double energy = 0;
    for(auto j : neigh_lists[i]){
        // get distance that is already computed
        rowvec dr = distance_dir_cube.tube(i, j);
        double dist_sqrd = distance_mat(i, j);

        // get types of particle j
        typej = box->types(j);
        double first_term = pow(dist_sqrd/sigma_mat(typei, typej), -3);
        double second_term = first_term * first_term;
        acc += 24 * epsilon_mat(typei, typej) * (2*second_term - first_term) * dr / dist_sqrd;
        if(comp_energy){
            energy += 4 * epsilon_mat(typei, typej) * (first_term - second_term);
        }
    }
    acc /= box->masses(i);

    return energy;
}


double LennardJones::eval_acc(const mat positions, mat& accs, const bool comp_energy)
{
    /* Evaluate force acting on all particles
     */
    // update distance matrix and neighbor lists
    ForceField::distance_matrix(positions);
    build_neigh_lists(positions);

    double energy = 0;
    accs = zeros(box->npar, box->ndim);
    for(int i=0; i<box->npar; i++){
        rowvec acc;
        energy += eval_acc_par(positions, i, acc, comp_energy);
        accs.row(i) = acc;
    }
    return energy;
}
