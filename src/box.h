#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <valarray>
#include <cassert>
#include <chrono>

//#include "tqdm/tqdm.h"

#include "io.h"
#include "dump.h"
#include "thermo.h"
//#include "init_velocity.h"
#include "rng/mersennetwister.h"
#include "forcefield/lennardjones.h"
//#include "integrator/velocityverlet.h"
#include "sampler/metropolis.h"
#include "boundary/stillinger.h"


class Box
{
public:
    Box(std::string working_dir_in="", double temp_in=1., double chempot_in=0.); 

    // methods
    void set_forcefield(class ForceField* forcefield_in);
    //void set_integrator(class Integrator* integrator_in);
    void set_sampler(class Sampler* sampler_in);
    void set_rng(class RandomNumberGenerator* rng_in);
    void set_boundary(class Boundary* boundary_in);
    //void set_velocity(class Velocity* velocity_in);
        
    void set_temp(const double temp_in);
    void set_chempot(const double chempot_in);
    void set_mass(const std::string label, const double mass);

    void add_move(class Moves* move, const double prob);
    void add_particle(class Particle *particle);
    void add_particle(const std::string label, const std::valarray<double> r);
    void add_particles(std::vector<class Particle *> particles);
    //void add_particles(const string chem_symbol, const mat positions_in);
    //void add_particles(const string chem_symbol, const mat positions_in, const mat velocities_in);
    //void add_particles(const vector<string> chem_symbols_in, const mat positions_in);
    //void add_particles(const vector<string> chem_symbols_in, const mat positions_in, const mat velocities_in);

    void snapshot(const std::string filename);
    void set_dump(const int freq, const std::string filename, const std::vector<std::string> outputs);
    void set_thermo(const int freq, const std::string filename, const std::vector<std::string> outputs);

    void check_particle_types();
    void init_simulation();
    int get_maxiter(const int nsteps);
    void print_info();

    void run_md(const int nsteps);
    void run_mc(const int nsteps, const int nmoves);

    // variables
    class ForceField* forcefield = nullptr;
    //class Integrator* integrator = nullptr;
    class Sampler* sampler = nullptr;
    class RandomNumberGenerator* rng = nullptr;
    class Dump* dump = nullptr;
    class Thermo* thermo = nullptr;
    class Boundary* boundary = nullptr;
    //class Velocity* velocity = nullptr;

    int npar, ndim, ntype, nmove, step, npar_prev;
    double temp, chempot, poteng, time;

    std::vector<class Particle *> particles;
    std::vector<std::string> unique_labels;
    std::vector<double> unique_masses;
    std::vector<class Moves*> moves;
    std::vector<double> moves_prob;
};
