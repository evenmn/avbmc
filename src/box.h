#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <valarray>

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
    Box(std::string  = "", double = 1., double = 0.); 

    // methods
    void set_forcefield(class ForceField*);
    //void set_integrator(class Integrator* integrator_in);
    void set_sampler(class Sampler*);
    void set_rng(class RandomNumberGenerator*);
    void set_boundary(class Boundary*);
    //void set_velocity(class Velocity* velocity_in);
        
    void set_temp(double);
    void set_chempot(double);
    void set_mass(std::string, double);

    void add_move(class Moves*, double);
    void add_particle(class Particle*);
    void add_particle(std::string, std::valarray<double>);
    void add_particles(std::vector<class Particle *>);
    void add_molecule_type(std::string, double);
    void add_molecule_type(std::vector<std::string>, double, double, std::vector<std::valarray<double> >, int = 0);

    void snapshot(std::string);
    void set_dump(int, std::string, std::vector<std::string>);
    void set_thermo(int, std::string, std::vector<std::string>);

    void check_particle_types();
    void init_simulation();
    int get_maxiter(int);
    void print_logo();
    void print_info();
    void print_mc_info();

    void run_md(int);
    void run_mc(int, int);

    // variables
    class ForceField* forcefield = nullptr;
    //class Integrator* integrator = nullptr;
    class Sampler* sampler = nullptr;
    class RandomNumberGenerator* rng = nullptr;
    class Dump* dump = nullptr;
    class Thermo* thermo = nullptr;
    class Boundary* boundary = nullptr;
    //class Velocity* velocity = nullptr;
    class MoleculeTypes* molecule_types = nullptr;

    int npar, ndim, ntype, nmove, nprocess, step;
    double temp, chempot, poteng, time;

    std::vector<class Particle *> particles;
    std::vector<std::string> unique_labels;
    std::vector<double> unique_masses;
    std::vector<class Moves*> moves;
    std::vector<double> moves_prob;
};
