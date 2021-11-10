#include "io.h"


vector<string> split(const string s)
{
    /* Split string by whitespace
     */
    stringstream ss(s);
    istream_iterator<string> begin(ss);
    istream_iterator<string> end;
    vector<string> vstrings(begin, end);
    return vstrings;
}

mat read_xyz(const string filename, vector<string>& chem_symbols)
{
    /* Read xyz-file consisting of one time frame
     */
    int npar, ndim;
    mat positions;
    ifstream f(filename);
    if(f.is_open()){
        string line;
        int line_num = 0;
        while(getline(f, line))
        {
            if(line_num == 0){
                // line containing number of particles
                npar = stoi(line);
            }
            else if(line_num == 1){
                // ignore info line
                line_num ++;
                continue;
            }
            else{ 
                // coordinate line
                vector<string> splitted = split(line);
                if(line_num == 2){
                    // count number of dimensions
                    ndim = splitted.size() - 1;
                    positions = zeros(npar, ndim);
                }
                chem_symbols.push_back(splitted[0]);
                for(int i=1; i<ndim; i++){
                    positions(line_num-2, i) = stod(splitted[i]);
                }
            }
            line_num ++;
        }
    }
    else{
        cout << "Could not open file '" + filename + "'! Aborting." << endl;
        exit(0);
    }
    return positions;
}

/*
void write_xyz(const string filename, const mat positions, const vector<string> chem_symbols, const string info, const bool append)
{
    // 
    //
    int npar = positions.n_rows;
    int ndim = positions.n_cols;

    ofstream f;
    if(append){
        f.open(filename, ofstream::out | ofstream::app);
    }
    else{
        f.open(filename, ofstream::out);
    }
    f << npar << endl;
    f << info << endl;;
    for(int i=0; i<npar; i++){
        f << chem_symbols[i];
        for(int j=0; j<ndim; j++){
            f << " " << positions(i, j);
        }
        f << endl;
    }
    f.close();
}
*/

void write_xyz(ofstream& f, const mat data, const vector<string> chem_symbols, const string info)
{
    /* 
     */
    int npar = data.n_rows;
    int ndim = data.n_cols;

    f << npar << endl;
    f << info << endl;;
    for(int i=0; i<npar; i++){
        f << chem_symbols[i];
        for(int j=0; j<ndim; j++){
            f << " " << data(i, j);
        }
        f << endl;
    }
}
