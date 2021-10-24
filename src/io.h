#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>


using namespace std;
using namespace arma;


mat read_xyz(const string filename);
void write_xyz(const string filename, const mat positions, const vector<string> chem_symbols, const string info = "", const bool append=false);
