#include <iostream>
#include <string>
#include <armadillo>

#include "core.h"
#include "initposition.h"
#include "initvelocity.h"
#include "sampler.h"
#include "moves.h"
#include "forcefield.h"

using namespace std;
using namespace arma;


int main()
{
    string working_directory = "simulation";
    vec position = {0., 0., 0.};
    sim Core(working_directory, position);
    sim.snapshot("initial.xyz");

    sim.dump(1, "dump.xyz", "x", "y", "z");

    // relax with molecular dynamics simulation
    sim.thermo(1, "md.log", "step", "time", "poteng", "kineng");
    sim.run_md(steps=100);
    sim.snapshot("after_md.xyz");

    // run Monte Carlo
    sim.set_sampler(Metropolis());
    sim.add_move(Teans(dx=0.01), 0.3);
    sim.add_move(TransMH(dc=0.01, Ddt=0.01), 0.7);
    sim.thermo(1, "mc.log", "step", "poteng", "acc_ratio");
    sim.run_mc(steps=1000, out="log");
    sim.snapshot("final.xyz");

    return 0;
}
