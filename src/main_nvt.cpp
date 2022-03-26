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
#include "init_position.h"
#include "rng/mersennetwister.h"
#include "boundary/open.h"
#include "boundary/periodic.h"
#include "forcefield/lennardjones.h"
//#include "sampler/umbrella.h"
#include "sampler/metropolis.h"
#include "moves/trans.h"
#include "moves/transmh.h"


int main()
{
    // initialize system
    System system("simulation");
    system.set_temp(0.7);
    MersenneTwister rng;
    system.set_rng(&rng);

    // initialize umbrella sampling with square function
    //double n_ = 32.;
    //double b = - 27/8 * std::pow(n_, 1/3.);
    //double a = 1000.0;
    //auto f = [a, b] (const int n) { return (a * (n + b * std::pow(n, 2/3.))); };
    //auto f = [] (const int n) { return (0.007 * (n - 32) * (n - 32)); };
    Metropolis sampler(&system);
    //Umbrella sampler(&system, f);
    system.set_sampler(&sampler);

    // initialize box
    Box box(&system, 2);
    system.add_box(&box);

    // set forcefield
    LennardJones forcefield(&box, "params.lj");
    box.set_forcefield(&forcefield);

    // set boundary
    Periodic boundary(&box, {30.0, 30.0, 30.0});
    //Open boundary(&box);
    box.set_boundary(&boundary);

    // add particles
    box.add_particles("Ar", fcc(6, 30));
    //box.add_particle("Ar", {0,0,0});
    //box.add_particle("Ar", {1.5,0,0});

    // initialize translation and AVBMC moves
    Trans move1(&system, &box, 0.1);
    TransMH move2(&system, &box, 0.1, 0.1);
    //system.add_move(&move1, 0.5);
    system.add_move(&move2, 0.5);

    // set sampling outputs
    box.set_dump(10, "mc.xyz", {"x", "y", "z"});
    box.set_thermo(100, "mc.log", {"step", "atoms", "poteng"});

    // run Monte Carlo simulation
    box.snapshot("initial.xyz");
    system.run_mc(3e3, 1);
    //box.snapshot("final.xyz");

    // dump number of status with a certain system size to file
    box.write_nsystemsize("nsystemsize.txt");
    return 0;
}
