#include <iostream>
#include <functional>
#include <memory>

#include "box.h"
#include "system.h"
#include "rng/mersennetwister.h"
#include "boundary/stillinger.h"
#include "forcefield/lennardjones.h"
#include "sampler/metropolis.h"
#include "sampler/umbrella.h"
#include "moves/moves.h"
#include "moves/trans.h"
//#include "moves/avbmc.h"
#include "moves/avbmcin.h"
#include "moves/avbmcout.h"
//#include "moves/avbmcinmol.h"


int main()
{
    // initialize system
    System system("simulation");
    system.set_temp(0.7);
    system.set_chempot(-1.3);
    LennardJones forcefield(&system, "params.lj");
    system.set_forcefield(&forcefield);
    MersenneTwister* rng = new MersenneTwister();
    system.set_rng(rng);

    // initialize umbrella sampling with square function
    //auto f = [] (const int n) { return (0.012 * (n - 32) * (n - 32)); };
    //Umbrella sampler(&system, f);
    Metropolis sampler(&system);
    system.set_sampler(&sampler);

    // initialize box with Stillinger boundary
    // criterion initialized with one Argon atom
    Box box(&system);
    Stillinger boundary(&box, 1.5);
    box.set_boundary(&boundary);
    box.add_particle("Ar", {0, 0, 0});
    //box.add_particle("Ar", {1.5, 0, 0});
    //box.add_particle("Ar", {0, 1.5, 0});
    //box.add_particle("Ar", {1.5, 1.5, 0});
    system.add_box(&box);

    // initialize translation and AVBMC moves
    Trans move1(&system, &box, 0.01);
    AVBMCIn move2(&system, &box, "Ar", 0.9, 1.5);
    system.add_move(&move1, 0.94);
    system.add_move(&move2, 0.06);
    //system.add_move(std::make_shared<AVBMCOut>(&system, box, 1.5), 0.5);

    // set sampling outputs
    //box.set_dump(1, "mc.xyz", {"x", "y", "z"});
    //box.set_thermo(1, "mc.log", {"step", "atoms", "poteng"});

    // run Monte Carlo simulation
    //box->snapshot("initial.xyz");
    //system.run_mc(10000, 1);
    //box.add_particle("Ar", {3.0, 0.0, 0.0});
    system.run_mc(10000, 1);
    //box.snapshot("final.xyz");

    // dump number of status with a certain system size to file
    //box.write_nsystemsize("nsystemsize.txt");

    delete rng;
    
    return 0;
}
