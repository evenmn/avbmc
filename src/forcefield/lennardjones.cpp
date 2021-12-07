#include "lennardjones.h"
#include "../box.h"
#include "../particle.h"
#include "../molecule.h"

LennardJones::LennardJones(Box* box_in)
    : ForceField(box_in)
{
    /* This is the default constructor
     * when parameter file is not given. It
     * only works for pure Argon
     */
    label1_vec = {"Ar"};
    label2_vec = {"Ar"};
    sigma_vec = {1.};
    epsilon_vec = {1.};
    rc_vec = {5.};
    nline = 1;
}

LennardJones::LennardJones(Box* box_in, const std::string params)
    : ForceField(box_in)
{
    read_param_file(params);
}

void LennardJones::read_param_file(const std::string params)
{
    /* Read parameter file and store parameters
     * globally. It takes the following form:
     * <label1> <label2> <sigma> <epsilon> <rc>
     */

    nline = 0;
    std::ifstream infile(params);
    std::string line, label1, label2;
    double sigma, epsilon, rc;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        if (line.rfind("#", 0) == 0){
            // comments are allowed in parameter file
            continue;
        }
        else if (iss >> label1 >> label2 >> sigma >> epsilon >> rc){ 
            label1_vec.push_back(label1);
            label2_vec.push_back(label2);
            sigma_vec.push_back(sigma);
            epsilon_vec.push_back(epsilon);
            rc_vec.push_back(rc);
            nline ++;
        }
        else{
            std::cout << "Warning: Corrupt line in parameter file" << std::endl;
            std::cout << "Line: " + line << std::endl;
        }
    }
}


void LennardJones::sort_params()
{
    /* Parameters have to be sorted with respect to
     * the particle types, and are stored in matrices.
     */

    // link list of chemical symbols to list of type indices
    std::vector<int> types1_vec;
    std::vector<int> types2_vec;
    for(std::string label : label1_vec){
        bool assigned = false;
        for(int j=0; j<box->ntype; j++){
            if(label == box->unique_labels[j]){
                types1_vec.push_back(j);
                assigned = true;
            }
        }
        assert(assigned);
    }
    for(std::string label : label2_vec){
        bool assigned = false;
        for(int j=0; j<box->ntype; j++){
            if(label == box->unique_labels[j]){
                types2_vec.push_back(j);
                assigned = true;
            }
        }
        assert(assigned);
    }

    // allocate memory for matrices
    sigma_mat = new double*[box->ntype];
    epsilon_mat = new double*[box->ntype];
    rc_sqrd_mat = new double*[box->ntype];
    shift_mat = new double*[box->ntype];
    for(int i=0; i<box->ntype; i++){
        sigma_mat[i] = new double[box->ntype];
        epsilon_mat[i] = new double[box->ntype];
        rc_sqrd_mat[i] = new double[box->ntype];
        shift_mat[i] = new double[box->ntype];
    }
    // fill up matrices with parameters
    int type1, type2;
    for(int i=0; i<nline; i++){
        type1 = types1_vec[i];
        type2 = types2_vec[i];
        sigma_mat[type1][type2] = sigma_vec[i];
        sigma_mat[type2][type1] = sigma_vec[i];
        epsilon_mat[type1][type2] = epsilon_vec[i];
        epsilon_mat[type2][type1] = epsilon_vec[i];
        double rc = rc_vec[i];
        rc_sqrd_mat[type1][type2] = rc * rc;
        rc_sqrd_mat[type2][type1] = rc * rc;
        shift_mat[type1][type2] = std::pow(rc, 12) - std::pow(rc, 6);
        shift_mat[type2][type1] = std::pow(rc, 12) - std::pow(rc, 6);
    }
}

/*
std::vector<int> LennardJones::build_neigh_list(const mat positions, const int i)
{
    // Build neighbor list for a particle i
    //

    // declare variables
    int typei, typej;
    vector<int> neigh_list;

    typei = box->particle_types[i];
    for(int j=0; j<i; j++){
        typej = box->particle_types[j];
        if(tmp_distance_mat(i, j) < rc_sqrd_mat(typei, typej)){
            neigh_list.push_back(j);
        }
    }
    for(int j=i+1; j<tmp_npar; j++){
        typej = box->particle_types[j];
        if(tmp_distance_mat(j, i) < rc_sqrd_mat(typej, typei)){
            neigh_list.push_back(j);
        }
    }
    return neigh_list;
}
*/
/*
void LennardJones::build_neigh_lists()
{
    // Build neigh lists for all particles. For Lennard
    // Jones, we will have one neighbor list only
    //
    for(int i=0; i<box->npar; i++){
        neigh_lists.push_back(build_neigh_list(box->positions, i));
    }
}
*/
/*
double LennardJones::eval_acc_element(const mat positions, const int i, const int j, rowvec& acc, const bool comp_energy)
{
    // Evaluate acceleration between two particles i and j
    //

    // declare variables
    int typei, typej;
    double first_term, second_term, dist_sqrd, energy;
    rowvec dr;

    // get distance that is already computed
    dr = distance_dir_cube.tube(i, j);
    dist_sqrd = distance_mat(i, j);


    // get types of particle j
    typei = box->particle_types[i];
    typej = box->particle_types[j];

    // calculate acceleration on particle i from particle j
    first_term = pow(dist_sqrd/sigma_mat(typei, typej), -3);
    second_term = first_term * first_term;
    acc += 24 * epsilon_mat(typei, typej) * (2*second_term - first_term) * dr / dist_sqrd / box->particle_masses[i];

    // calculate interactio energy between the two particles
    if(comp_energy){
        energy = 4 * epsilon_mat(typei, typej) * (first_term - second_term);
    }
    return energy;
}
*/
/*
double LennardJones::eval_acc_par(const mat positions, const int i, rowvec& acc, const bool comp_energy)
{
    // Evaluate force acting on particle i. This is needed only when
    // particle i is moved, and therefore we also have to update the
    // distance between particle i and all other particles.
    //

    // update neighbor list of particle i
    build_neigh_list(positions, i);
    
    acc.zeros(box->ndim);
    double energy = 0;
    for(int j : neigh_lists[i]){
        energy += eval_acc_element(positions, i, j, acc, comp_energy);
    }
    return energy;
}
*/
/*
std::vector<class Particle *> LennardJones::build_neigh_list(const std::vector<class Particle *> particles, const int i)
{
    // declare variables
    int typei = particles[i].type; 
    double sdist;
    std::valarray<double> posi = particles[i].r;
    std::vector<class Particle *> neigh_list;

    
    for(Particle particle : particles){
        sdist = std::sum(std::pow(particle.r - ref_pos, 2));
        if(sdist < rc_sqrd_mat[typei][particle.type]){
            neigh_list.push_back(particle);
        }
    }
    return neigh_list;
}
*/


double LennardJones::comp_energy_mol(const std::vector<Particle *> particles, Molecule* molecule)
{
    /* Compute energy contribution from a molecule
     */
    double energy = 0.;
    for(int i : molecule->atoms_idx){
        energy += comp_energy_par(particles, i);
    }
    return energy;
}


double LennardJones::comp_energy_par(const std::vector<Particle *> particles, const int i)
{
    /*
     */
    // declare variables
    int typei = particles[i]->type; 
    double sdist, ss6, ss12;

    double energy = 0.;
    for(int j=0; j<i; j++){
        sdist = std::pow(particles[j]->r - particles[i]->r, 2).sum();
        int typej = particles[j]->type;
        if(sdist < rc_sqrd_mat[typei][typej]){
            ss6 = std::pow(sigma_mat[typei][typej] / sdist, 3);
            ss12 = ss6 * ss6;
            energy += epsilon_mat[typei][typej] * (ss12 - ss6 - shift_mat[typei][typej]);
        }
    }
    for(int j=i+1; j<box->npar; j++){
        sdist = std::pow(particles[j]->r - particles[i]->r, 2).sum();
        int typej = particles[j]->type;
        if(sdist < rc_sqrd_mat[typei][typej]){
            ss6 = std::pow(sigma_mat[typei][typej] / sdist, 3);
            ss12 = ss6 * ss6;
            energy += epsilon_mat[typei][typej] * (ss12 - ss6 - shift_mat[typei][typej]);
        }
    }
    return (4 * energy);
}

/*
double LennardJones::update_force_par(const mat positions, const int i)
{
    // Update force matrix and potential energy matrix for a particle.
    // This is needed by Monte Carlo simulations.
    //
    
    // declare
    int typei, typej;
    double dist_sqrd_inv, s6, s12;
    rowvec dr;

    cout << "update_force_par1" << endl;
    std::vector<int> tmp_neigh_list = build_neigh_list(positions, i);

    cout << "update_force_par2" << endl;
    typei = tmp_particle_types[i];
    // compute force
    for(int j : tmp_neigh_list){
        dr = tmp_distance_dir_cube.tube(j, i);
        cout << "update_force_par3" << endl;
        dist_sqrd_inv = (1.0) / tmp_distance_mat(j, i);
        cout << "update_force_par4" << endl;
        typej = tmp_particle_types[j];
        cout << "update_force_par5" << endl;
        s6 = pow(sigma_mat(typei, typej) * dist_sqrd_inv, 3);
        cout << "update_force_par6" << endl;
        s12 = s6 * s6;
        cout << "update_force_par7" << endl;
        tmp_force_mat.tube(j, i) = 24 * epsilon_mat(typei, typej) * (2*s6 - s12) * dr * dist_sqrd_inv;
        cout << "update_force_par8" << endl;
        tmp_poteng_mat(j, i) = 4 * epsilon_mat(typei, typej) * (s12 - s6);
        cout << "update_force_par9" << endl;
    }

    // update acceleration
    mat forces_par = sum(force_mat, 1);
    tmp_accelerations.copy_size(forces_par);
    for(int i = 0; i<tmp_npar; i++){
        tmp_accelerations.row(i) = forces_par.row(i) / tmp_particle_masses[i];
    }

    return sum(sum(tmp_poteng_mat));
}
*/
/*
double LennardJones::update_force_all()
{
    // Update force on all particles and the energy contribution of
    // each particle pair. Called every time for molecular dynamics
    // simulations.
    //

    // declare variables
    int typei, typej;
    double s6, s12, dist_sqrd_inv;
    rowvec dr;

    // update distance matrix and neighbor lists
    update_distance_all();
    build_neigh_lists();

    // compute force
    for(int i=0; i<box->npar; i++){
        typei = box->particle_types[i];
        for(int j=0; j<i; j++){
            dr = distance_dir_cube.tube(j, i);
            dist_sqrd_inv = (1.0) / distance_mat(j, i);
            typej = box->particle_types[j];
            s6 = pow(sigma_mat(typei, typej) * dist_sqrd_inv, 3);
            s12 = s6 * s6;
            force_mat.tube(j, i) = 24 * epsilon_mat(typei, typej) * (2*s6 - s12) * dr * dist_sqrd_inv;
            poteng_mat(j, i) = 4 * epsilon_mat(typei, typej) * (s12 - s6);
        }
    }

    // update acceleration
    mat forces_par = sum(force_mat, 1);
    for(int i = 0; i<forces_par.n_rows; i++){
        box->accelerations.row(i) = forces_par.row(i) / box->particle_masses[i];
    }

    return sum(sum(poteng_mat));
}
*/
/*
double LennardJones::eval_acc(const mat positions, mat& accs, vec& potengs, const bool comp_energy)
{
    // Evaluate forces acting on all particles
    //
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
*/
/*
double LennardJones::comp_force_par(const rowvec pos, rowvec &acc)
{
    // Evaluate the energy change of the system
    // when a particle is added at position pos
    //
    double energy_cum = 0;
    acc.zeros(box->ndim);
    for(int i=0; i<box->npar; i++){
        rowvec dr = box->positions.row(i) - pos;
        double dist2 = as_scalar(dr * dr.as_col());
        double dist_inv6 = pow(dist2, -3);           // / sigma
        double dist_inv12 = pow(dist_inv6, 2);
        energy_cum += 4 * (dist_inv12 - dist_inv6);  // * epsilon
        acc += 24 * (2 * dist_inv6 - dist_inv12) * dr / dist2;
    }
    return energy_cum;
}
*/

LennardJones::~LennardJones()
{
    for(int i = 0; i < box->ntype; i++){
        delete (sigma_mat[i]);
        delete (epsilon_mat[i]);
        delete (rc_sqrd_mat[i]);
        delete (shift_mat[i]);
    }
    delete (sigma_mat);
    delete (epsilon_mat);
    delete (rc_sqrd_mat);
    delete (shift_mat);
}

