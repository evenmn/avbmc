#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <armadillo>

#include "box.h"

using namespace std;
using namespace arma;


void parser(int argc, char** argv[]);
void set(Box& box, string keyword);
void add(Box& box, string keyword);
void take(Box& box, string keyword);
void run(Box& box, string keyword);
