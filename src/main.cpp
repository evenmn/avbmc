#include <iostream>
#include <string>
#include <armadillo>

#include "box.h"
//#include "forcefield/forcefield.h"
#include "forcefield/lennardjones.h"
//#include "integrator/integrator.h"
#include "integrator/euler.h"

using namespace std;
using namespace arma;


int main()
{
    // initialize box
    string working_directory = "simulation";
    double temp = 300.;
    double chempot = 0.;
    Box box(working_directory, temp, chempot);

    // initialize particles
    int type = 0;
    double mass = 1.;
    mat position = randn(10, 3);
    mat velocity = randn(10, 3);
    string chem = "Ar";
    box.add_particles(type, mass, position, velocity, chem);

    box.snapshot("initial.xyz");

    //box.dump(1, "dump.xyz", "x", "y", "z");

    // relax with molecular dynamics simulation
    box.set_integrator(new Euler(&box));
    box.set_forcefield(new LennardJones(&box, ".in"));
    //box.thermo(1, "md.log", "step", "time", "poteng", "kineng");
    box.run_md(100);
    box.snapshot("after_md.xyz");
    /*
    // run Monte Carlo
    box.set_sampler(Metropolis());
    box.add_move(Teans(dx=0.01), 0.3);
    box.add_move(TransMH(dc=0.01, Ddt=0.01), 0.7);
    box.thermo(1, "mc.log", "step", "poteng", "acc_ratio");
    box.run_mc(steps=1000, out="log");
    box.snapshot("final.xyz");
    */
    return 0;
}
