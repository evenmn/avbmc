#include <iostream>
#include <string>
#include <valarray>
#include <vector>

#include "box.h"
#include "init_position.h"
#include "forcefield/lennardjones.h"
//#include "integrator/euler.h"
#include "sampler/metropolis.h"
#include "moves/trans.h"
//#include "moves/transmh.h"
#include "moves/avbmc.h"

using namespace std;


int main()
{
    // initialize box
    Box box("simulation", 300., 0.);

    // initialize particles
    std::valarray<double> r1(1, 3);
    box.add_particle("Ar", r1);
    //std::valarray<double> r2(2, 3);
    //box.add_particle("Ar", r2);
    box.set_mass("Ar", 1.);

    box.snapshot("initial.xyz");
    
    // run Monte Carlo
    box.set_sampler(new Metropolis(&box));
    box.add_move(new Trans(&box, 0.1), 0.94);
    box.add_move(new AVBMC(&box, 0.9, 1.5), 0.06);
    std::vector<std::string> outputs_dump = {"xyz"}; //, "acceptance_ratio"};
    box.set_dump(1, "mc.xyz", outputs_dump);
    std::vector<std::string> outputs_thermo = {"step", "atoms", "poteng"}; //, "acceptance_ratio"};
    box.set_thermo(1, "mc.log", outputs_thermo);
    box.run_mc(100000, 1);
    //box.snapshot("final.xyz");
    
    return 0;
}
