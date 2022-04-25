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

#include "boundary/open.h"
#include "box.h"
#include "forcefield/idealgas.h"
#include "forcefield/vashishta.h"
#include "io.h"
#include "moves/trans.h"
#include "particle.h"
#include "rng/mersennetwister.h"
#include "sampler/metropolis.h"
#include "sampler/umbrella.h"
#include "system.h"
//#include "moves/avbmcmolin.h"
#include "moves/avbmcmolinres.h"
//#include "moves/avbmcmolout.h"
#include "moves/avbmcmoloutres.h"
#include "constraint/minneigh.h"
#include "constraint/maxneigh.h"
#include "constraint/mindistance.h"
#include "constraint/maxdistance.h"
#include "constraint/stillinger.h"


int main()
{
    // ======  initialize system  ======
    System system("simulation");
    system.set_temp(320.);
    system.set_chempot(10.);
    MersenneTwister rng;
    system.set_rng(&rng);

    // ======  initialize sampler  ======
    auto f = [] (const int n) { return (0.003 * (n - 20) * (n - 20)); };
    //Umbrella sampler(&system, f);
    Metropolis sampler(&system);
    system.set_sampler(&sampler);

    // ======  initialize box  ======
    // initialize box with open boundaries
    Box box(&system, 2);
    system.add_box(&box);
    Open boundary(&box);
    box.set_boundary(&boundary);

    // ===== initialize forcefield =====
    Vashishta forcefield(&box, "H2O.nordhagen.vashishta");
    // IdealGas forcefield(&box, {"H", "O"});
    box.set_forcefield(&forcefield);

    // ======  initialize box constraints  ======
    Stillinger constraint(&box);
    constraint.set_criterion("O", "O", 3.0);
    constraint.set_criterion("O", "H", 1.6);
    constraint.set_criterion("H", "H", 0.0);
    MinNeigh constraint1(&box, "O", "O", 3.0, 2);
    // MinNeigh constraint2(&box, "O", "H", 1.3, 2);
    // MaxNeigh constraint3(&box, "O", "H", 1.0, 2);
    // MaxDistance constraint1(&box, "O", "O", 4.0);
    // MinDistance constraint2(&box, "O", "O", 3.0);
    // MaxDistance constraint2(&box, "O", "H", 1.6);
    // MaxDistance constraint3(&box, "H", "H", 0.0);
    box.add_constraint(&constraint);
    box.add_constraint(&constraint1);
    //box.add_constraint(&constraint2);
    //box.add_constraint(&constraint3);

    // ======  initialize particle positions  ======
    std::vector<Particle> single_molecule = read_xyz("water.xyz");
    //std::vector<Particle> bulk_water = read_xyz("water_bulk.xyz");
    box.add_particles(single_molecule);

    // ======  initialize moves  ======
    Trans move1(&system, &box, 0.1);
    //Trans move2(&system, &box, 1.0);
    AVBMCMolInRes move3(&system, &box, single_molecule, 3., 0.9, 4.0);
    AVBMCMolOutRes move4(&system, &box, single_molecule, 4.0, 3.);
    system.add_move(&move1, 0.5);
    //system.add_move(&move2, 0.45);
    system.add_move(&move3, 0.25);
    system.add_move(&move4, 0.25);

    // set sampling outputs
    box.set_dump(1, "mc.xyz", {"x", "y", "z"});
    box.set_thermo(1, "mc.log", {"step", "atoms", "poteng"});

    // run Monte Carlo simulation
    //box.snapshot("initial.xyz");
    system.run_mc(1e4, 1);
    box.snapshot("final.xyz");
    //system.run_mc(1e5, 1);

    // dump number of status with a certain system size to file
    box.write_nsystemsize("nsystemsize_biased_water.txt");
    return 0;
}
