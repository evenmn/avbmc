#include <iostream>
#include <string>
#include <vector>
#include <armadillo>

#include "forcefield.h"
#include "moves.h"
#include "sampler.h"
#include "dump.h"
#include "thermo.h"

using namespace std;
using namespace arma;

class Core{
    public:
        Core(string working_dir_in, double temp_in, double chempot_in);

        set_forcefield(ForceField forcefield_in);
        set_integrator(Integrator integrator_in);
        set_sampler(Sampler sampler_in);
        
        add_move(Moves move);
        add_particles(int type, double mass, mat position, mat velocity);

        snapshot(string filename);
        dump(int freq, string filename, vector<string> outputs);
        thermo(int freq, string filename, vector<string> outputs);

    privat:
        write_xyz(string filename);
        void Core::write_xyz(string filename, mat arma_mat, string chem="X", string info="type x y z", bool overwrite=false)

        ForceField forcefield = nullptr;
        Integrator integrator = nullptr;
        Sampler sampler = nullptr;

        double temp, chempot;

        vector<int> types;
        vector<double> masses;
        vector<mat> positions;
        vector<mat> velocities;
        vector<string> chem_symbols;
        vector<Moves> moves;
        vector<double> moves_prob;
};
