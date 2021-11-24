#include <iostream>
#include <string>
#include <valarray>
#include <vector>

#include "box.h"
#include "forcefield/lennardjones.h"
#include "sampler/umbrella.h"
#include "moves/trans.h"
#include "moves/avbmc.h"

using namespace std;


int main()
{
    // initialize box
    Box box("simulation");
    box.set_temp(0.7);
    box.set_chempot(-1.3);

    // initialize particle at (0, 0, 0)
    box.add_particle("Ar", {0, 0, 0});
    box.set_mass("Ar", 1.);
    box.set_forcefield(new LennardJones(&box, "params.lj"));

    // initialize Monte Carlo
    double k = 0.007;
    double nc = 32.;
    auto f = [nc, k] (const int n) { return (k * (n - nc) * (n - nc)); };
    box.set_sampler(new Umbrella(&box, f));
    box.add_move(new Trans(&box, 0.01), 0.94);
    box.add_move(new AVBMC(&box, 0.9, 1.5), 0.06);

    // set outputs
    box.set_dump(1, "mc.xyz", {"x", "y", "z"});
    box.set_thermo(1, "mc.log", {"step", "atoms", "poteng", "acceptanceratio"});

    // run Monte Carlo simulation
    box.snapshot("initial.xyz");
    box.run_mc(10000, 1);
    box.snapshot("final.xyz");
    
    return 0;
}
