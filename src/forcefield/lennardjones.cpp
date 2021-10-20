#include "lennardjones.h"

LennardJones::LennardJones(string params)
{
    read_param_file(params);
}

void LennardJones::read_param_file(string params)
{
    /* Read parameter file and store parameters
     * globally. It takes the following form:
     * <type1> <type2> <sigma> <epsilon> <rc>
     */
    sigma = 1;
    epsilon = 1;
    rc = 3;
    rc_sqrd = rc * rc;
}


void LennardJones::build_neigh_list(const mat positions, const int i)
{
    /* Build neighbor list for a particle i
     */
    vector<int> neigh_list;
    for(int j=0; j<i; j++){
        if(distance_mat(i, j) < rc_sqrd){
            neigh_list.push_back(j);
        }
    }
    neigh_lists(i) = neigh_list;
}


void LennardJones::build_neigh_lists(const mat positions)
{
    /* Build neigh lists for all particles
     */
    neigh_lists(npar);
    for(int i=0; i<npar; i++){
        build_neigh_list(positions, i);
    }
}


double LennardJones::eval_force_par(const mat positions, const int i, vec& force, const bool comp_energy=false)
{
    /* Evaluate force acting on particle i. This is needed only when
     * particle i is moved, and therefore we also have to update the
     * distance between particle i and all other particles.
     */
    ForceField::distance_matrix_cross(positions, i);
    build_neigh_list(positions, i);
    
    force(ndim, fill:zeros);
    double energy = 0;
    for(auto j : neigh_lists(i)){
        vec dr = distance_dir_cube.tube(i, j);
        double dist_sqrd = distance_matrix(i, j);
        double first_term = pow(dist_sqrd/sigma, -3);
        double second_term = first_term * first_term;
        force += (2*second_term - first_term) * dr / dist_sqrd;
        if(comp_energy){
            energy += first_term - second_term
        }
    }
    force *= 24 * epsilon;
    energy *= 4 * epsilon;
    return energy;
}


double LennardJones::eval_force(const mat positions, mat& forces, const bool comp_energy=false)
{
    /* Evaluate force acting on all particles
     */
    double energy = 0;
    forces(npar, ndim);
    for(int i=0; i<npar; i++){
        vec force;
        energy += eval_force_particle(positions, i, force, comp_energy);
        forces.row(i) = force;
    }
    return energy;
}
