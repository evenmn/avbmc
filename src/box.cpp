#include <iostream>
#include <string>
#include <vector>
#include <valarray>
#include <cassert>
#include <chrono>

#include "box.h"
#include "dump.h"
#include "moves/moves.h"
#include "particle.h"
#include "molecule.h"


/* -------------------------------------------------------
   Box constructor 
---------------------------------------------------------- */

Box::Box(std::string working_dir_in, double temp_in, double chempot_in)
{
    //working_dir = working_dir_in;
    temp = temp_in;
    chempot = chempot_in;

    time = poteng = 0.;
    npar = ntype = nmove = nprocess = step = 0;
    ndim = 3;

    rng = new MersenneTwister();
    //integrator = new VelocityVerlet(this);
    forcefield = new LennardJones(this);
    sampler = new Metropolis(this);
    boundary = new Stillinger(this);
    //velocity = new Zero();
    molecule_types = new MoleculeTypes(this);

    std::vector<std::string> outputs;
    dump = new Dump(this, 0, "", outputs);
    outputs = {"step", "atoms", "poteng"};
    thermo = new Thermo(this, 1, "", outputs);
}


/* --------------------------------------------------------
   Set box temperature
----------------------------------------------------------- */

void Box::set_temp(const double temp_in)
{
    temp = temp_in;
}


/* --------------------------------------------------------
   Set chemical potential of system
----------------------------------------------------------- */

void Box::set_chempot(const double chempot_in)
{
    chempot = chempot_in;
}


/* --------------------------------------------------------
   Set mass of chemical symbol. Masses of all chemical symbols
   have to be given, as the software does not look up the
   masses in a table.
----------------------------------------------------------- */

void Box::set_mass(const std::string label, const double mass)
{
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


/* ------------------------------------------------------
   Add new molecule type when molecule is defined by
   one atom.
--------------------------------------------------------- */

void Box::add_molecule_type(std::string element, const double molecule_prob)
{
    std::vector<std::string> elements = {element};
    std::valarray<double> default_atom(0.0, ndim);
    std::vector<std::valarray<double> > default_mol = {default_atom};
    molecule_types->add_molecule_type(elements, 0.0, 0, molecule_prob, default_mol);
}


/* ------------------------------------------------------
   Add new molecule type consisting of several atoms
--------------------------------------------------------- */

void Box::add_molecule_type(std::vector<std::string> elements, const double rc, const double molecule_prob,
                            std::vector<std::valarray<double> > default_mol, const int com_atom)
{
    molecule_types->add_molecule_type(elements, rc, com_atom, molecule_prob, default_mol);
}


/* ------------------------------------------------------
   Dump snapshot of system using the "write_xyz"-function. 
--------------------------------------------------------- */
   
void Box::snapshot(const std::string filename)
{
    std::vector<std::string> outputs = {"xyz"};
    Dump* tmp_dump = new Dump(this, 1, filename, outputs);
    tmp_dump->print_frame();
    delete tmp_dump;
}

/* -----------------------------------------------------
   Specify dump output
-------------------------------------------------------- */

void Box::set_dump(const int freq, const std::string filename, const std::vector<std::string> outputs)
{
    dump = new Dump(this, freq, filename, outputs);
}


/* -----------------------------------------------------
   Specify thermo output
-------------------------------------------------------- */

void Box::set_thermo(const int freq, const std::string filename, const std::vector<std::string> outputs)
{
    thermo = new Thermo(this, freq, filename, outputs);
}


/* -----------------------------------------------------
   The particles in the system have to be a subset of 
   the particles that are given mass and  type. That 
   has to be checked after parsing, but before 
   simulation is started.
-------------------------------------------------------- */

void Box::check_particle_types()
{
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

    // Sort forcefield parameters according to particle types
    forcefield->sort_params();

    // If molecule configuration is not set, let all single particles
    // be considered as molecules
    if (!molecule_types->configured)
    {
        for(std::string label : unique_labels){
            add_molecule_type(label, 1./ntype);
        }
    }

    // Since particles are associated with types rather than labels,
    // also molecule configuration labels have to be converted to types
    for(std::vector<std::string> elements : molecule_types->molecule_elements){
        std::vector<int> types;
        for (int i=0; i < ntype; i++){
            for(std::string element : elements){
                if(element == unique_labels[i]){
                    types.push_back(i);
                }
            }
        }
        molecule_types->molecule_types.push_back(types); 
    }
}


/* -------------------------------------------------------
   Initialize variables needed before simulation
---------------------------------------------------------- */
void Box::init_simulation()
{
    thermo->print_header();
}


/* -------------------------------------------------------
   Returns the last iteration
---------------------------------------------------------- */

int Box::get_maxiter(const int nsteps)
{
    int maxiter;
    if(step == 0){
        maxiter = nsteps + 1;
    }
    else{
        maxiter = nsteps + step;
    }
    return maxiter;
}


/* -------------------------------------------------------
   Print logo header
---------------------------------------------------------- */

void Box::print_logo()
{
    std::cout << std::endl;
    std::cout << " █████╗ ██╗   ██╗██████╗ ███╗   ███╗ ██████╗" << std::endl;
    std::cout << "██╔══██╗██║   ██║██╔══██╗████╗ ████║██╔════╝" << std::endl;
    std::cout << "███████║██║   ██║██████╔╝██╔████╔██║██║     " << std::endl;
    std::cout << "██╔══██║╚██╗ ██╔╝██╔══██╗██║╚██╔╝██║██║     " << std::endl;
    std::cout << "██║  ██║ ╚████╔╝ ██████╔╝██║ ╚═╝ ██║╚██████╗" << std::endl;
    std::cout << "╚═╝  ╚═╝  ╚═══╝  ╚═════╝ ╚═╝     ╚═╝ ╚═════╝" << std::endl;
    std::cout << std::endl;
}


/* -------------------------------------------------------
   Print information about simulation
---------------------------------------------------------- */

void Box::print_info()
{
    std::time_t start_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::cout << "Started computation at " << std::ctime(&start_time)
              << "Running on " << nprocess << " CPU threads using MPI" << std::endl;

    //std::cout << std::fixed;
    //std::cout << std::boolalpha;
    //std::cout << std::setprecision(6);
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "            Box Information " << std::endl;
    std::cout << "=========================================" << std::endl;
    std::cout << "Number of particles:      " << npar << std::endl;
    std::cout << "Number of dimensions:     " << ndim << std::endl;
    std::cout << "Number of particle types: " << ntype << std::endl;
    //std::cout << "Forcefield:               " << forcefield->label << std::endl;
    //std::cout << "Random number generator:  " << rng->label << std::endl;
    std::cout << std::endl;
    std::cout << "Unique particle types:";
    for(int i=0; i < ntype; i++){
        std::cout << " " << unique_labels[i];
    }
    std::cout << std::endl;
    for(int i=0; i < molecule_types->ntype; i++){
        std::cout << "  Molecule type " + std::to_string(i+1) + ":" << std::endl;
        std::cout << "    Atoms:            ";
        for(std::string element : molecule_types->molecule_elements[i]){
            std::cout << " " << element;
        }
        std::cout << std::endl;
        std::cout << "    Critical distance: " << std::to_string(molecule_types->rcs[i]) << std::endl;
        std::cout << "    Probability:       " << std::to_string(molecule_types->molecule_probs[i]) << std::endl;
    }
    std::cout << std::endl;
}


/* --------------------------------------------------------
   Print Monte Carlo information
----------------------------------------------------------- */

void Box::print_mc_info()
{
    std::cout << std::endl;
    std::cout << "           Monte Carlo information " << std::endl;
    std::cout << "=========================================" << std::endl;
    std::cout << "Temperature:              " << temp << std::endl;
    std::cout << "Chemical potential:       " << chempot << std::endl;
    //std::cout << "Sampling method:          " << sampler->label << std::endl;
    std::cout << "Number of move types:     " << nmove << std::endl;
    std::cout << std::endl;
    for(int i=0; i < nmove; i++){
        std::cout << "  Move type " + std::to_string(i+1) + ": ";
        std::cout << moves[i]->repr();
        std::cout << "    Probability: " << std::to_string(moves_prob[i]) << std::endl;
    }
    std::cout << std::endl;
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

/* -------------------------------------------------------
   Run Monte Carlo simulation
---------------------------------------------------------- */

void Box::run_mc(const int nsteps, const int nmoves)
{
    print_logo();
    print_info();
    print_mc_info();
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
