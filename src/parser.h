#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <functional> 
#include <cassert>
#include <iterator>
#include <vector>

#include "io.h"
#include "box.h"
#include "init_position.h"
//#include "init_velocity.h"

//#include "integrator/integrator.h"
//#include "integrator/euler.h"
//#include "integrator/eulercromer.h"
//#include "integrator/velocityverlet.h"
//#include "integrator/rungekutta4.h"

#include "forcefield/forcefield.h"
#include "forcefield/lennardjones.h"
#include "forcefield/vashishta.h"

#include "sampler/sampler.h"
#include "sampler/metropolis.h"
#include "sampler/umbrella.h"

#include "moves/moves.h"
#include "moves/trans.h"
//#include "moves/transmh.h"
#include "moves/avbmc.h"

#include "boundary/boundary.h"
//#include "boundary/fixed.h"
#include "boundary/stillinger.h"


using namespace std;


void parser(int argc, char** argv);
void set(Box& box, const std::vector<std::string> splitted, const int argc);
void add(Box& box, const std::vector<std::string> splitted, const int argc);
void take(Box& box, const std::vector<std::string> splitted, const int argc);
void run(Box& box, const std::vector<std::string> splitted, const int argc);
void thermo(Box& box, const std::vector<std::string> splitted, const int argc);
void dump(Box& box, const std::vector<std::string> splitted, const int argc);
