/* ----------------------------------------------------------------------------
   General info, authors, url, license etc...
------------------------------------------------------------------------------- */


/* ----------------------------------------------------------------------------
   This file contains the main function that is compiled when building with
   'make dev' or 'make debug'. It may be edited, but please make sure that you
   know what you do. The main policy is that all the objects used in the code
   are created here, and their *adresses* are passed to the various functions.
   Let's take the Metropolis sampler object as an example:

    Method 1: 
    Metropolis sampler(&system);
    system.set_sampler(&sampler);

    Method 2:
    Metropolis* sampler = new Metropolis(&system);
    system.set_sampler(sampler);
    ...
    delete sampler;

    Avoid:
    system.set_sampler(new Metropolis(&system));

   The last example should be avoided, as it will lead to an unpleasant memory
   leak (the pointer is never deleted).
------------------------------------------------------------------------------- */

#include <iostream>
#include <functional>

#include "box.h"
#include "system.h"
#include "rng/mersennetwister.h"
#include "boundary/stillinger.h"
#include "forcefield/lennardjones.h"
#include "sampler/metropolis.h"
#include "moves/trans.h"
#include "moves/avbmcin.h"
#include "moves/avbmcout.h"


int main()
{
    // initialize system
    System system("simulation");
    system.set_temp(0.741);
    system.set_chempot(-0.31512);
    LennardJones forcefield(&system, "params.lj");
    system.set_forcefield(&forcefield);
    MersenneTwister rng;
    system.set_rng(&rng);

    // initialize umbrella sampling with square function
    Metropolis sampler(&system);
    system.set_sampler(&sampler);

    // initialize box with Stillinger boundary
    // criterion initialized with one Argon atom
    Box box(&system);
    system.add_box(&box);
    Stillinger boundary(&box, 1.5);
    box.set_boundary(&boundary);
    box.add_particle("Ar", {0, 0, 0});

    // initialize translation and AVBMC moves
    Trans move1(&system, &box, 0.1);
    AVBMCIn move2(&system, &box, "Ar", 0.9, 1.5);
    AVBMCOut move3(&system, &box, "Ar", 1.5);
    system.add_move(&move1, 0.94);
    system.add_move(&move2, 0.03);
    system.add_move(&move3, 0.03);

    // set sampling outputs
    //box.set_dump(1, "mc.xyz", {"x", "y", "z"});
    //box.set_thermo(1, "mc.log", {"step", "atoms", "poteng"});

    // run Monte Carlo simulation
    box.snapshot("initial.xyz");
    system.run_mc(1e7, 1);
    box.snapshot("final.xyz");

    // dump number of status with a certain system size to file
    box.write_nsystemsize("nsystemsize.txt");
    return 0;
}
