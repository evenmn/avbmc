#include <iostream>
#include <string>
#include <armadillo>
#include <vector>

#include "box.h"
#include "init_position.h"
#include "forcefield/lennardjones.h"
#include "integrator/euler.h"
#include "sampler/metropolis.h"
#include "moves/trans.h"
#include "moves/transmh.h"
#include "moves/avbmcout.h"

using namespace std;
using namespace arma;


int main()
{
    // initialize box
    Box box("simulation", 300., 0.);

    // initialize particles
    int type = 0;
    double mass = 1.;
    mat position = fcc(3, 5., 3);
    mat velocity = randn(position.n_rows, position.n_cols);
    string chem = "Ar";
    box.add_particles(type, mass, position, velocity, chem);

    box.snapshot("initial.xyz");

    //box.dump(1, "dump.xyz", "x", "y", "z");

    // relax with molecular dynamics simulation
    box.set_integrator(new Euler(&box));
    box.set_forcefield(new LennardJones(&box, ".in"));
    vector<string> outputs = {"Step", "PotEng", "AcceptanceRatio"};
    box.set_thermo(2, "md.log", outputs);
    outputs = {"xyz"};
    box.set_dump(100, "mc.xyz", outputs);
    //box.run_md(100);
    //box.snapshot("after_md.xyz");
    
    // run Monte Carlo
    box.set_sampler(new Metropolis(&box));
    box.add_move(new Trans(&box, 0.1), 0.5);
    box.add_move(new TransMH(&box, 0.01, 0.01), 0.5);
    //box.add_move(new AVBMCOut(&box, 0.5, 2.), 0.25);
    //box.thermo(1, "mc.log", "step", "poteng", "acc_ratio");
    //box.run_mc(steps=1000, out="log");
    box.run_mc(1000, 100);
    box.snapshot("final.xyz");
    
    return 0;
}
