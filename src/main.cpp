#include <iostream>
#include <string>
#include <valarray>
#include <vector>

#include "box.h"
#include "system.h"
#include "forcefield/lennardjones.h"
#include "sampler/umbrella.h"
#include "moves/trans.h"
#include "moves/avbmc.h"


int main()
{
    // initialize system
    System system("simulation");
    system.set_temp(0.7);
    system.set_chempot(-1.3);
    system.set_forcefield(new LennardJones(&box, "params.lj"));

    // initialize Monte Carlo
    double k = 0.007;
    double nc = 32.;
    auto f = [nc, k] (const int n) { return (k * (n - nc) * (n - nc)); };
    system.set_sampler(new Umbrella(&box, f));
    system.add_move(new Trans(&box, 0.01), 0.94);
    system.add_move(new AVBMC(&box, 0.9, 1.5), 0.06);

    // initialize box
    Box box1(system);
    box1.add_particle("Ar", {0, 0, 0});

    // set outputs
    box1.set_dump(1, "mc.xyz", {"x", "y", "z"});
    box1.set_thermo(1, "mc.log", {"step", "atoms", "poteng", "acceptanceratio"});

    // run Monte Carlo simulation
    box1.snapshot("initial.xyz");
    system.run_mc(10000, 1);
    box1.snapshot("final.xyz");
    
    return 0;
}
