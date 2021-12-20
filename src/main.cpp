#include <iostream>
#include <string>
#include <valarray>
#include <vector>

#include "box.h"
#include "system.h"
#include "boundary/stillinger.h"
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
    system.set_forcefield(new LennardJones(&system, "params.lj"));

    // initialize Monte Carlo and umbrella sampling
    double k = 0.007;
    double nc = 32.;
    auto f = [nc, k] (const int n) { return (k * (n - nc) * (n - nc)); };
    system.set_sampler(new Umbrella(&system, f));

    // initialize box
    Box box1(&system);
    box1.set_boundary(new Stillinger(&box1, 1.5));
    box1.add_particle("Ar", {0, 0, 0});

    // initialize moves
    system.add_move(new Trans(&system, &box1, 0.01), 1.);
    system.add_move(new AVBMC(&system, &box1, 0.9, 1.5), 0.06);

    // set outputs
    box1.set_dump(1, "mc.xyz", {"x", "y", "z"});
    box1.set_thermo(1, "mc.log", {"step", "atoms", "poteng"});

    // run Monte Carlo simulation
    box1.snapshot("initial.xyz");
    //system.run_mc(10000, 1);
    box1.snapshot("final.xyz");
    
    return 0;
}
