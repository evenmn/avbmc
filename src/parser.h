#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cassert>
#include <armadillo>
#include <iterator>
#include <vector>

#include "io.h"
#include "box.h"
#include "init_position.h"
#include "init_velocity.h"

#include "integrator/integrator.h"
#include "integrator/euler.h"
#include "integrator/eulercromer.h"
#include "integrator/velocityverlet.h"
#include "integrator/rungekutta4.h"

#include "forcefield/forcefield.h"
#include "forcefield/lennardjones.h"

#include "sampler/sampler.h"
#include "sampler/metropolis.h"

#include "moves/moves.h"
#include "moves/trans.h"
#include "moves/transmh.h"
#include "moves/avbmcout.h"

#include "boundary/boundary.h"
#include "boundary/fixed.h"
#include "boundary/stillinger.h"


using namespace std;
using namespace arma;


void parser(int argc, char** argv);
void set(Box& box, const vector<string> splitted, const int argc);
void add(Box& box, const vector<string> splitted, const int argc);
void take(Box& box, const vector<string> splitted, const int argc);
void run(Box& box, const vector<string> splitted, const int argc);
void thermo(Box& box, const vector<string> splitted, const int argc);
void dump(Box& box, const vector<string> splitted, const int argc);
//vector<string> split(const string s);
