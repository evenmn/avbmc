#include "core.h"

Core::Core(string working_dir_in, double temp_in, double chempot_in){
    working_dir = working_dir_in;
    temp = temp_in;
    chempot = chempot_in;
}

void Core::set_forcefield(ForceField forcefield_in)
{
    forcefield = forcefield_in;
}

void Core::set_integrator(Integrator integrator_in)
{
    integrator = integrator_in;
}

void Core::set_sampler(Sampler sampler_in)
{
    sampler = sampler_in;
}

void Core::add_move(Moves move, double prob)
{
    moves.push_back(move);
    moves_prob.push_back(prob);
}

void Core::add_particles(int type, double mass, mat position, mat velocity, string chem)
{
    /* Add particles of type "type", mass "mass",
     * positions "position" and velocities "velocity"
     */
    types.push_back(type);
    masses.push_back(mass);
    position.push_back(position);
    velocities.push_back(velocity);
    chem_symbols.push_back(chem);
}

void Core::write_xyz(string filename, mat arma_mat, string chem, string info, bool overwrite)
{
    /* Write an Armadillo matrix containing
     * positions at a given timestep to file with
     * xyz-format. The cube has to have particles on
     * first axis and dimensions on second axis
     */
    if(overwrite){
        ofstream f;
    }
    else{
        ifstream f;
    }
    f.open(filename);
    int npar = arma_cube.n_rows;
    int ndim = arma_cube.n_cols;
    f << npar << endl;
    f << info << endl;
    for(int i=0; i<npar; i++){
        f << chem;
        for(int j=0; j<ndim; j++){
            f << " " << arma_mat(i, j) << endl;
        }
    }
    f.close();
}

void Core::snapshot(string filename){
    /* Dump snapshot of system using the
     * "write_xyz"-function. 
     */
    write_xyz(filename, positions, chem_symbols, "type x y z", true);
}
