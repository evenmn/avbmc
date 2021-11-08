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
     * <chemsymbol1> <chemsymbol2> <sigma> <epsilon> <rc>
     */

    chem_symbol1_vec = {"Ar"};
    chem_symbol2_vec = {"Ar"};
    sigma_vec = {1.};
    epsilon_vec = {1.};
    rc_vec = {3.};
    nline = 1;

    /*
    ifstream infile(params);
    string line;
    int type1, type2;
    double sigma, epsilon, rc;
    while (getline(infile, line))
    {
        istringstream iss(line);
        if (line[0] == "#"){
            continue;
        }
        else if (!(iss >> type1 >> type2 >> sigma >> epsilon >> rc)){ 
            break; 
        }
    }
    */
    /*
    int type1, type2;
    double sigma, epsilon, rc;
    while (infile >> type1 >> type2 >> sigma >> epsilon >> rc)
    {
        type1_vec.push_back(type1);
        type2_vec.push_back(type2);
        sigma_vec.push_back(sigma);
        epsilon_vec.push_back(epsilon);
        rc_vec.push_back(rc);
    }
    */
}

void LennardJones::sort_params()
{
    /* Parameters have to be sorted with respect to
     * the particle types, and are stored in matrices.
     */

    // link list of chemical symbols to list of type indices
    vector<int> types1_vec;
    vector<int> types2_vec;
    for(string chem_symbol : chem_symbol1_vec){
        bool assigned = false;
        for(int j=0; j<box->ntype; j++){
            if(chem_symbol == box->unique_chem_symbols[j]){
                types1_vec.push_back(j);
                assigned = true;
            }
        }
        assert(assigned);
    }
    for(string chem_symbol : chem_symbol2_vec){
        bool assigned = false;
        for(int j=0; j<box->ntype; j++){
            if(chem_symbol == box->unique_chem_symbols[j]){
                types1_vec.push_back(j);
                assigned = true;
            }
        }
        assert(assigned);
    }

    // fill up matrices with parameters
    sigma_mat.zeros(box->ntype, box->ntype);
    epsilon_mat.zeros(box->ntype, box->ntype);
    rc_sqrd_mat.zeros(box->ntype, box->ntype);
    
    int type1, type2;
    for(int i=0; i<nline; i++){
        type1 = types1_vec[i];
        type2 = types2_vec[i];
        sigma_mat(type1, type2) = sigma_vec[i];
        sigma_mat(type2, type1) = sigma_vec[i];
        epsilon_mat(type1, type2) = epsilon_vec[i];
        epsilon_mat(type2, type1) = epsilon_vec[i];
        double rc = rc_vec[i];
        rc_sqrd_mat(type1, type2) = rc * rc;
        rc_sqrd_mat(type2, type1) = rc * rc;
    }
}


vector<int> LennardJones::build_neigh_list(const mat positions, const int i)
{
    /* Build neighbor list for a particle i
     */

    // declare variables
    int typei, typej;
    vector<int> neigh_list;

    // update distances between particle i and all other particles
    ForceField::distance_matrix_cross(positions, i);

    typei = box->types(i);
    for(int j=0; j<i; j++){
        typej = box->types(j);
        if(distance_mat(i, j) < rc_sqrd_mat(typei, typej)){
            neigh_list.push_back(j);
        }
    }
    for(int j=i+1; j<box->npar; j++){
        typej = box->types(j);
        if(distance_mat(j, i) < rc_sqrd_mat(typej, typei)){
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


double LennardJones::eval_acc_element(const mat positions, const int i, const int j, rowvec& acc, const bool comp_energy)
{
    /* Evaluate acceleration between two particles i and j
     */

    // declare variables
    int typei, typej;
    double first_term, second_term, dist_sqrd, energy;
    rowvec dr;

    // get distance that is already computed
    dr = distance_dir_cube.tube(i, j);
    dist_sqrd = distance_mat(i, j);


    // get types of particle j
    typei = box->types(i);
    typej = box->types(j);

    // calculate acceleration on particle i from particle j
    first_term = pow(dist_sqrd/sigma_mat(typei, typej), -3);
    second_term = first_term * first_term;
    acc += 24 * epsilon_mat(typei, typej) * (2*second_term - first_term) * dr / dist_sqrd / box->masses(i);

    // calculate interactio energy between the two particles
    if(comp_energy){
        energy = 4 * epsilon_mat(typei, typej) * (first_term - second_term);
    }
    return energy;
}


double LennardJones::eval_acc_par(const mat positions, const int i, rowvec& acc, const bool comp_energy)
{
    /* Evaluate force acting on particle i. This is needed only when
     * particle i is moved, and therefore we also have to update the
     * distance between particle i and all other particles.
     */

    // update neighbor list of particle i
    build_neigh_list(positions, i);
    
    acc.zeros(box->ndim);
    double energy = 0;
    for(int j : neigh_lists[i]){
        energy += eval_acc_element(positions, i, j, acc, comp_energy);
    }
    return energy;
}


double LennardJones::eval_acc(const mat positions, mat& accs, vec& potengs, const bool comp_energy)
{
    /* Evaluate forces acting on all particles
     */
    double energy_cum, energy_par;
    energy_cum = 0;
    potengs.zeros(box->npar);
    accs.zeros(box->npar, box->ndim);
    for(int i=0; i<box->npar; i++){
        rowvec acc;
        energy_par = eval_acc_par(positions, i, acc, comp_energy);
        potengs(i) = energy_par;
        energy_cum += energy_par;
        accs.row(i) = acc;
    }
    return energy_cum;
}
