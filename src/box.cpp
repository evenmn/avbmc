#include "box.h"
#include "particle.h"

Box::Box(string working_dir_in, double temp_in, double chempot_in)
{
    //working_dir = working_dir_in;
    temp = temp_in;
    chempot = chempot_in;

    time = 0.;
    npar = ntype = nmove = step = 0;

    rng = new MersenneTwister();
    //integrator = new VelocityVerlet(this);
    forcefield = new LennardJones(this);
    sampler = new Metropolis(this);
    boundary = new Stillinger(this);
    //velocity = new Zero();

    std::vector<std::string> outputs;
    dump = new Dump(this, 0, "", outputs);
    outputs = {"step", "atoms", "poteng"};
    thermo = new Thermo(this, 1, "", outputs);
}

void Box::set_temp(const double temp_in)
{
    /* Set system temperature. 
     */
    temp = temp_in;
}

void Box::set_chempot(const double chempot_in)
{
    /* Set chemical potential of system.
     */
    chempot = chempot_in;
}

void Box::set_mass(const std::string label, const double mass)
{
    /* Set mass of chemical symbol. Masses of all chemical symbols
     * have to be given, as the software does not look up the
     * masses in a table.
     */
    ntype ++;
    unique_labels.push_back(label);
    unique_masses.push_back(mass);
}

void Box::set_forcefield(class ForceField* forcefield_in)
{
    forcefield = forcefield_in;
}
/*
void Box::set_integrator(class Integrator* integrator_in)
{
    integrator = integrator_in;
}
*/
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
/*
void Box::set_velocity(class Velocity* velocity_in)
{
    velocity = velocity_in;
    velocities = velocity->get_velocity(npar, ndim);
}
*/
void Box::add_move(class Moves* move, double prob)
{
    /* Add move type and the corresponding probability.
     * The sum of the probabilities have to be 1.
     */
    nmove ++;
    moves.push_back(move);
    moves_prob.push_back(prob);
}

void Box::add_particle(class Particle* particle)
{
    /* Add single particle
     */
    npar ++;
    ndim = particle->r.size();
    particles.push_back(particle);
}

void Box::add_particle(const std::string label, const std::valarray<double> r)
{
    /*
     */
    npar ++;
    ndim = r.size();
    Particle *particle = new Particle(label, r);
    particles.push_back(particle);
}

void Box::add_particles(std::vector<Particle *> particles_in)
{
    /*
     */
    npar += particles_in.size();
    ndim = particles_in[0]->r.size();
    particles.insert(particles.end(), particles_in.begin(), particles_in.end());
}

/*
void Box::add_particles(const string chem_symbol, const mat positions_in)
{
    // Add particles of the same type, chemical symbol 'chem_symbol' at positions
    // 'positions_in'. No initial velocities.
    //

    // check dimensionality
    int npar_added = positions_in.n_rows;
    npar += npar_added;
    ndim = positions_in.n_cols;
    
    // join old and new particles positions and velocities
    positions = join_cols(positions, positions_in);
    velocities = join_cols(velocities, zeros(npar_added, ndim));

    // join old and new particles chemical symbols
    vector<string> chem_symbol_vec(npar_added, chem_symbol);
    chem_symbols.insert(chem_symbols.end(), chem_symbol_vec.begin(), chem_symbol_vec.end());
}

void Box::add_particles(const string chem_symbol, const mat positions_in, const mat velocities_in)
{
    // Add particles of the same type, chemical symbol 'chem_symbol' at positions
    // 'positions_in' and velocities 'velocities_in'.
    //

    // check dimensionality
    int npar_added = positions_in.n_rows;
    npar += npar_added;
    ndim = positions_in.n_cols;
    assert (velocities_in.n_rows == npar_added);
    
    // join old and new particles positions and velocities
    positions = join_cols(positions, positions_in);
    velocities = join_cols(velocities, velocities_in);

    // join old and new particles chemical symbols
    vector<string> chem_symbol_vec(npar_added, chem_symbol);
    chem_symbols.insert(chem_symbols.end(), chem_symbol_vec.begin(), chem_symbol_vec.end());
}

void Box::add_particles(const vector<string> chem_symbols_in, const mat positions_in)
{
    // Add particles of potentially various chemical symbols at positions
    // 'positions_in'. No initial velocities.
    //

    // check dimensionality
    int npar_added = positions_in.n_rows;
    npar += npar_added;
    ndim = positions_in.n_cols;
    
    // join old and new particles positions and velocities
    positions = join_cols(positions, positions_in);
    velocities = join_cols(velocities, zeros(npar_added, ndim));

    // join old and new particles chemical symbols
    chem_symbols.insert(chem_symbols.end(), chem_symbols_in.begin(), chem_symbols_in.end());
}

void Box::add_particles(const vector<string> chem_symbols_in, const mat positions_in, const mat velocities_in)
{
    // Add particles of potentially various chemical symbols at positions
    // 'positions_in' and velocities 'velocities_in'.
    //

    // check dimensionality
    int npar_added = positions_in.n_rows;
    npar += npar_added;
    ndim = positions_in.n_cols;
    assert (velocities_in.n_rows == npar_added);
    
    // join old and new particles positions and velocities
    positions = join_cols(positions, positions_in);
    velocities = join_cols(velocities, velocities_in);

    // join old and new particles chemical symbols
    chem_symbols.insert(chem_symbols.end(), chem_symbols_in.begin(), chem_symbols_in.end());
}
*/

void Box::snapshot(const std::string filename){
    /* Dump snapshot of system using the
     * "write_xyz"-function. 
     */

    std::vector<std::string> outputs = {"xyz"};
    class Dump* tmp_dump = new Dump(this, 1, filename, outputs);
    tmp_dump->print_frame();
    delete tmp_dump;
}


void Box::set_dump(const int freq, const std::string filename, const std::vector<std::string> outputs)
{
    /* Specify dump output
     */

    dump = new Dump(this, freq, filename, outputs);
}


void Box::set_thermo(const int freq, const std::string filename, const std::vector<std::string> outputs)
{
    /* Specify thermo output
     */

    thermo = new Thermo(this, freq, filename, outputs);
}


void Box::check_particle_types()
{
    /* The particles in the system have to be a subset
     * of the particles that are given mass and 
     * type. That has to be checked after parsing,
     * but before simulation is started.
     */

    // Check that all particles are assigned a mass
    bool not_assigned_mass_all = 0;
    for(Particle *particle : particles){
        bool assigned_mass = false;
        for(int j=0; j<ntype; j++){
            if(particle->label == unique_labels[j]){
                particle->type = j;
                particle->mass = unique_masses[j];
                assigned_mass = true;
            }
        }
        not_assigned_mass_all += !assigned_mass;
    }
    assert(!not_assigned_mass_all);

    // Check that all particles are assigned parameters
    forcefield->sort_params();
}


void Box::init_simulation()
{
    /* Initialize variables needed before simulation.
     */
    thermo->print_header();
}

int Box::get_maxiter(const int nsteps)
{
    /* Returns the last iteration
     */
    int maxiter;
    if(step == 0){
        maxiter = nsteps + 1;
    }
    else{
        maxiter = nsteps + step;
    }
    return maxiter;
}

void Box::print_info()
{
    /* Print information about simulation
     */
    cout << "Num particles: " << npar << endl;
    cout << "Num dimensions: " << ndim << endl;
    cout << "Num moves: " << nmove << endl;
    cout << "Num types: " << ntype << endl;
}

/*
void Box::run_md(const int nsteps)
{
    // Run molecular dynamics simulation
    //
    print_info();
    check_particle_types();
    init_simulation();

    int maxiter = get_maxiter(nsteps);

    // run molecular dynamics simulation
    auto t1 = chrono::high_resolution_clock::now();
    while(step<maxiter){
        time = step * integrator->dt;
        dump->print_frame();
        thermo->print_line();
        integrator->next_step();
        step ++;
    }
    auto t2 = chrono::high_resolution_clock::now();
    double duration_seconds = chrono::duration<double>(t2 - t1).count();
    cout << "Elapsed time: " << duration_seconds << "s" << endl;
}
*/

void Box::run_mc(const int nsteps, const int nmoves)
{
    /* Run Monte Carlo simulation
     */
    print_info();
    check_particle_types();
    init_simulation();

    int maxiter = get_maxiter(nsteps);

    // run Monte Carlo simulation
    auto t1 = chrono::high_resolution_clock::now();
    while(step < maxiter){
        dump->print_frame();
        thermo->print_line();
        sampler->sample(nmoves);
        step ++;
    }
    auto t2 = chrono::high_resolution_clock::now();
    double duration_seconds = chrono::duration<double>(t2 - t1).count();
    cout << "Elapsed time: " << duration_seconds << "s" << endl;
}
