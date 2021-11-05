#include "box.h"

Box::Box(string working_dir_in, double temp_in, double chempot_in)
{
    working_dir = working_dir_in;
    temp = temp_in;
    chempot = chempot_in;

    npar = 0;

    integrator = new Euler(this);
    forcefield = new LennardJones(this, ".in");
    sampler = new Metropolis(this);
    boundary = new Stillinger(this);
    rng = new MersenneTwister();
}

void Box::set_temp(double temp_in)
{
    temp = temp_in;
}

void Box::set_chempot(double chempot_in)
{
    chempot = chempot_in;
}

void Box::set_forcefield(class ForceField* forcefield_in)
{
    forcefield = forcefield_in;
}

void Box::set_integrator(class Integrator* integrator_in)
{
    integrator = integrator_in;
}


void Box::set_sampler(class Sampler* sampler_in)
{
    sampler = sampler_in;
}

void Box::set_rng(class RandomNumberGenerator* rng_in)
{
    rng = rng_in;
}

void Box::set_boundary(class Boundary* boundary_in)
{
    boundary = boundary_in;
}


void Box::add_move(class Moves* move, double prob)
{
    moves.push_back(move);
    moves_prob.push_back(prob);
}


void Box::add_particles(const int type, const double mass, const mat position, const mat velocity, const string chem)
{
    /* Add particles of type "type", mass "mass",
     * positions "position" and velocities "velocity"
     */

    // check dimensionality
    int npar_added = position.n_rows;
    npar += npar_added;
    ndim = position.n_cols;
    assert (velocity.n_rows == npar_added);
    
    // join old and new particles
    types = join_cols(types, type * ones(npar_added));
    masses = join_cols(masses, mass * ones(npar_added));
    positions = join_cols(positions, position);
    velocities = join_cols(velocities, velocity);

    vector<string> chem_symbol_vec(npar_added, chem);
    chem_symbol.insert(chem_symbol.end(), chem_symbol_vec.begin(), chem_symbol_vec.end());
}


void Box::snapshot(string filename){
    /* Dump snapshot of system using the
     * "write_xyz"-function. 
     */

    vector<string> outputs = {"xyz"};
    class Dump* tmp_dump = new Dump(this, 1, filename, outputs);
    tmp_dump->print_frame();
    delete tmp_dump;
}


void Box::set_dump(int freq, const string filename, const vector<string> outputs)
{
    /* Specify dump output
     */

    dump = new Dump(this, freq, filename, outputs);
}


void Box::set_thermo(int freq, const string filename, const vector<string> outputs)
{
    /* Specify thermo output
     */

    thermo = new Thermo(this, freq, filename, outputs);
}


void Box::run_md(int nsteps)
{
    /* Run molecular dynamics simulation
     */

    // compute initial acceleration (this should be done when adding particles?)
    forcefield->distance_mat = zeros(npar, npar);
    forcefield->distance_dir_cube = zeros(npar, npar, ndim);
    forcefield->build_neigh_lists(positions);
    poteng = forcefield->eval_acc(positions, accelerations, potengs, true);

    // run molecular dynamics simulation
    for(step=0; step<nsteps; step++){
        cout << step << endl;
        integrator->next_step();
        //time = step * integrator->dt;
        // dump
        cout << step << endl;
        write_xyz("dump.xyz", positions, chem_symbol, "", true);
    }
}


void Box::run_mc(int nsteps, int nmoves)
{
    /* Run Monte Carlo simulation
     */

    // compute initial acceleration (this should be done when adding particles?)
    forcefield->distance_mat = zeros(npar, npar);
    forcefield->distance_dir_cube = zeros(npar, npar, ndim);
    forcefield->build_neigh_lists(positions);
    poteng = forcefield->eval_acc(positions, accelerations, potengs, true);
    thermo->print_header();

    // run Monte Carlo simulation
    auto t1 = chrono::high_resolution_clock::now();
    for(step=0; step<nsteps; step++){
        //cout << step << endl;
        dump->print_frame();
        thermo->print_line();
        sampler->sample(nmoves);
        //write_xyz("dump.xyz", positions, chem_symbol, "", true);
    }
    auto t2 = chrono::high_resolution_clock::now();
    double duration_seconds = chrono::duration<double>(t2 - t1).count();
    cout << duration_seconds << endl;
}
