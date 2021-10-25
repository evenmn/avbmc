#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <armadillo>
#include <cassert>

//#include "tqdm/tqdm.h"

#include "io.h"
#include "rng/mersennetwister.h"
#include "forcefield/forcefield.h"
#include "forcefield/lennardjones.h"
#include "integrator/integrator.h"
#include "integrator/euler.h"
#include "moves/moves.h"
#include "sampler/sampler.h"
//#include "dump.h"
//#include "thermo.h"

using namespace std;
using namespace arma;

class Box
{
public:
    Box(string working_dir_in, double temp_in, double chempot_in); 

    void set_forcefield(class ForceField* forcefield_in);
    void set_integrator(class Integrator* integrator_in);
    void set_sampler(class Sampler* sampler_in);
    void set_rng(class RandomNumberGenerator* rng_in);
        
    void add_move(class Moves* move, const double prob);
    void add_particles(const int type, const double mass, const mat position, const mat velocity, const string chem);

    void snapshot(const string filename);
    //void dump(int freq, string filename, vector<string> outputs);
    //void thermo(int freq, string filename, vector<string> outputs);

    void run_md(int nsteps);
    void run_mc(int nsteps, int nmoves);

    string working_dir;

    class ForceField* forcefield = nullptr;
    class Integrator* integrator = nullptr;
    class Sampler* sampler = nullptr;
    class RandomNumberGenerator* rng = nullptr;

    int npar, ndim;
    double temp, chempot, poteng;

    vec types, masses;
    mat positions, velocities, accelerations;

    vector<string> chem_symbol;
    vector<class Moves*> moves;
    vector<double> moves_prob;
};
