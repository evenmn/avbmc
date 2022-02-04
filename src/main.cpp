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

#include "io.h"
#include "box.h"
#include "system.h"
#include "particle.h"
#include "rng/mersennetwister.h"
#include "boundary/stillinger.h"
#include "forcefield/vashishta.h"
#include "sampler/umbrella.h"
#include "sampler/metropolis.h"
#include "moves/trans.h"
#include "moves/avbmcmolin.h"
#include "moves/avbmcmolinres.h"
#include "moves/avbmcmolout.h"
#include "moves/avbmcmoloutres.h"


int main()
{
    // initialize system
    System system("simulation");
    system.set_temp(320.);
    system.set_chempot(-2.1);
    Vashishta forcefield(&system, "H2O.nordhagen.vashishta");
    system.set_forcefield(&forcefield);
    MersenneTwister rng;
    system.set_rng(&rng);

    // initialize umbrella sampling with square function
    auto f = [] (const int n) { return (0.003 * (n - 20) * (n - 20)); };
    Umbrella sampler(&system, f);
    //Metropolis sampler(&system);
    system.set_sampler(&sampler);

    // initialize box with Stillinger boundary
    // criterion initialized with one Argon atom
    Box box(&system);
    system.add_box(&box);
    Stillinger boundary(&box);
    boundary.set_crit("O", "O", 4.0);
    boundary.set_crit("O", "H", 1.6);
    boundary.set_crit("H", "H", 0.0);
    std::vector<Particle> single_molecule = read_xyz("water.xyz");
    //std::vector<Particle> bulk_water = read_xyz("water_bulk.xyz");
    box.set_boundary(&boundary);
    box.add_particles(single_molecule);

    // initialize translation and AVBMC moves
    Trans move1(&system, &box, 0.3);
    AVBMCMolInRes move2(&system, &box, single_molecule, 1.2, 0.9, 4.0);
    AVBMCMolOutRes move3(&system, &box, single_molecule, 4.0, 1.2);
    system.add_move(&move1, 0.94);
    system.add_move(&move2, 0.03);
    system.add_move(&move3, 0.03);

    // set sampling outputs
    box.set_dump(100, "mc.xyz", {"x", "y", "z"});
    box.set_thermo(100, "mc.log", {"step", "atoms", "poteng"});

    // run Monte Carlo simulation
    //box.snapshot("initial.xyz");
    system.run_mc(2e6, 1);
    box.snapshot("final.xyz");

    // dump number of status with a certain system size to file
    box.write_nsystemsize("nsystemsize_biased_water.txt");
    return 0;
}
