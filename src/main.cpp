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
#include "forcefield/vashishta.h"
#include "sampler/umbrella.h"
#include "moves/moves.h"
#include "moves/trans.h"
#include "moves/transmh.h"
#include "moves/avbmcin.h"
#include "moves/avbmcout.h"
#include "moves/avbmcmolin.h"


int main()
{
    // initialize system
    System system("simulation");
    system.set_temp(0.7);
    system.set_chempot(-1.3);
    //LennardJones forcefield(&system, "params.lj");
    Vashishta forcefield(&system, "SiC.vashishta");
    system.set_forcefield(&forcefield);
    MersenneTwister* rng = new MersenneTwister();
    system.set_rng(rng);

    // initialize umbrella sampling with square function
    auto f = [] (const int n) { return (0.012 * (n - 32) * (n - 32)); };
    Umbrella sampler(&system, f);
    system.set_sampler(&sampler);

    // initialize box with Stillinger boundary
    // criterion initialized with one Argon atom
    Box box(&system);
    Stillinger boundary(&box, 4.0);
    box.set_boundary(&boundary);
    box.add_particle("Si", {0, 0, 0});
    box.add_particle("C", {3, 0, 0});
    system.add_box(&box);

    // initialize translation and AVBMC moves
    Trans move1(&system, &box, 1.0);
    //TransMH move2(&system, &box, 1.0, 0.);
    //AVBMCIn move3(&system, &box, "Si", 0.9, 4.0);
    AVBMCMolIn move3(&system, &box, 0.9, 4.0);
    AVBMCOut move4(&system, &box, 4.0);
    system.add_move(&move1, 0.94);
    //system.add_move(&move2, 0.47);
    system.add_move(&move3, 0.03);
    system.add_move(&move4, 0.03);

    // set sampling outputs
    //box.set_dump(1, "mc.xyz", {"x", "y", "z"});
    //box.set_thermo(1, "mc.log", {"step", "atoms", "poteng"});

    // run Monte Carlo simulation
    //box->snapshot("initial.xyz");
    system.run_mc(100000, 1);
    //box.snapshot("final.xyz");

    // dump number of status with a certain system size to file
    box.write_nsystemsize("nsystemsize.txt");

    delete rng;
    return 0;
}
