#include <iostream>
#include <string>
#include <valarray>
#include <vector>

#include "box.h"
#include "forcefield/lennardjones.h"
#include "sampler/umbrella.h"
#include "moves/trans.h"
#include "moves/avbmcmol.h"


int main()
{
    // initialize box
    Box box("simulation");
    box.set_temp(0.7);
    box.set_chempot(-2.2);

    // initialize particle at (0, 0, 0)
    box.add_particle("Ar", {0, 0, 0});
    box.set_forcefield(new LennardJones(&box, "params.lj"));

    // initialize Monte Carlo
    double k = 0.012;
    double nc = 32.;
    auto f = [nc, k] (const int n) { return (k * (n - nc) * (n - nc)); };
    box.set_sampler(new Umbrella(&box, f));
    box.add_move(new Trans(&box, 0.1), 0.66);
    box.add_move(new AVBMCMol(&box, 0.9, 1.5), 0.34);

    // set outputs
    box.set_dump(6000, "mc.xyz", {"x", "y", "z"});
    box.set_thermo(1, "mc.log", {"step", "atoms", "poteng", "acceptanceratio"});

    // run Monte Carlo simulation
    box.run_mc(600000, 100);
    return 0;
}
