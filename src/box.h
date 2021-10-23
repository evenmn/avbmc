#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <armadillo>
#include <cassert>

#include "forcefield/forcefield.h"
#include "forcefield/lennardjones.h"
#include "integrator/integrator.h"
#include "integrator/euler.h"
//#include "moves/moves.h"
//#include "sampler/sampler.h"
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
    //void set_sampler(Sampler sampler_in);
        
    //void add_move(Moves move);
    void add_particles(const int type, const double mass, const mat position, const mat velocity, const string chem);

    void snapshot(string filename);
    void dump(int freq, string filename, vector<string> outputs);
    void thermo(int freq, string filename, vector<string> outputs);

    void run_md(int steps);

    string working_dir;

    class ForceField* forcefield = nullptr;
    class Integrator* integrator = nullptr;
    //Sampler sampler;

    int npar, ndim;
    double temp, chempot;

    vec types, masses;
    mat positions, velocities, accelerations;

    vector<string> chem_symbol;
    //vector<Moves> moves;
    //vector<double> moves_prob;

private:
    void write_xyz(string filename);
    //void Core::write_xyz(string filename, mat arma_mat, string chem="X", string info="type x y z", bool overwrite=false)

};
