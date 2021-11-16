#include <iostream>
#include <string>
#include <vector>

#include "box.h"
#include "init_position.h"
#include "forcefield/lennardjones.h"
//#include "integrator/euler.h"
#include "sampler/metropolis.h"
#include "moves/trans.h"
//#include "moves/transmh.h"
#include "moves/avbmcout.h"

using namespace std;


int main()
{
    // initialize box
    Box box("simulation", 300., 0.);

    // initialize particles
    Particle *particle;
    particle->label = "Ar";
    particle->r = {0., 0., 0.};
    box.add_particle(particle);
    box.set_mass("Ar", 1.);

    box.snapshot("initial.xyz");
    
    // run Monte Carlo
    box.set_sampler(new Metropolis(&box));
    box.add_move(new Trans(&box, 0.1), 1.0);
    //box.add_move(new AVBMCOut(&box, 0.5, 2.), 0.25);
    //box.thermo(1, "mc.log", "step", "poteng", "acc_ratio");
    box.run_mc(1000, 100);
    //box.snapshot("final.xyz");
    
    return 0;
}
