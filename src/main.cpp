#include <iostream>
#include <functional>

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

    // initialize umbrella sampling with square function
    auto f = [] (const int n) { return (0.007 * (n - 32) * (n - 32)); };
    system.set_sampler(new Umbrella(&system, f));

    // initialize box with Stillinger boundary
    // criterion initialized with one Argon atom
    Box box(&system);
    box.set_boundary(new Stillinger(&box, 1.5));
    box.add_particle("Ar", {0, 0, 0});
    system.add_box(&box);

    // initialize translation and AVBMC moves
    system.add_move(new Trans(&system, &box, 0.01), 0.94);
    system.add_move(new AVBMC(&system, &box, 0.9, 1.5), 0.06);

    // set sampling outputs
    //box.set_dump(1, "mc.xyz", {"x", "y", "z"});
    //box.set_thermo(1, "mc.log", {"step", "atoms", "poteng"});

    // run Monte Carlo simulation
    //box.snapshot("initial.xyz");
    system.run_mc(1000000, 1);
    //box.snapshot("final.xyz");

    // dump number of status with a certain system size to file
    box.write_nsystemsize("nsystemsize.txt");
    
    return 0;
}
