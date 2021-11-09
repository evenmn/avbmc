#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <iterator>
#include <armadillo>


using namespace std;
using namespace arma;

vector<string> split(const string s);
mat read_xyz(const string filename, vector<string>& chem_symbols);
//void write_xyz(const string filename, const mat positions, const vector<string> chem_symbols, const string info = "", const bool append=false);
void write_xyz(ofstream& f, const mat data, const vector<string> chem_symbols, const string info);

