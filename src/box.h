#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <armadillo>
#include <cassert>
#include <chrono>

//#include "tqdm/tqdm.h"

#include "io.h"
#include "dump.h"
#include "thermo.h"
#include "rng/mersennetwister.h"
//#include "forcefield/forcefield.h"
#include "forcefield/lennardjones.h"
//#include "integrator/integrator.h"
#include "integrator/euler.h"
//#include "moves/moves.h"
//#include "sampler/sampler.h"
#include "sampler/metropolis.h"
#include "boundary/stillinger.h"

using namespace std;
using namespace arma;


class Box
{
public:
    Box(string working_dir_in="", double temp_in=1., double chempot_in=0.); 

    void set_forcefield(class ForceField* forcefield_in);
    void set_integrator(class Integrator* integrator_in);
    void set_sampler(class Sampler* sampler_in);
    void set_rng(class RandomNumberGenerator* rng_in);
    void set_boundary(class Boundary* boundary_in);
        
    void set_temp(const double temp_in);
    void set_chempot(const double chempot_in);
    void set_mass(const string chem, const double mass);

    void add_move(class Moves* move, const double prob);
    void add_particles(const string chem_symbol, const mat positions_in);
    void add_particles(const string chem_symbol, const mat positions_in, const mat velocities_in);
    void add_particles(const vector<string> chem_symbols_in, const mat positions_in);
    void add_particles(const vector<string> chem_symbols_in, const mat positions_in, const mat velocities_in);

    void snapshot(const string filename);
    void set_dump(int freq, string filename, vector<string> outputs);
    void set_thermo(int freq, string filename, vector<string> outputs);

    void check_particle_types();
    void init_simulation();

    void run_md(int nsteps);
    void run_mc(int nsteps, int nmoves);

    string working_dir;

    class ForceField* forcefield = nullptr;
    class Integrator* integrator = nullptr;
    class Sampler* sampler = nullptr;
    class RandomNumberGenerator* rng = nullptr;
    class Dump* dump = nullptr;
    class Thermo* thermo = nullptr;
    class Boundary* boundary = nullptr;

    int npar, ndim, ntype, step;
    double temp, chempot, poteng, time;

    vec types, masses, potengs;
    mat positions, velocities, accelerations;

    vector<string> chem_symbols;
    vector<int> particle_types;
    vector<string> unique_chem_symbols;
    vector<double> unique_masses;
    vector<class Moves*> moves;
    vector<double> moves_prob;
};
