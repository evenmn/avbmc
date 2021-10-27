#include "io.h"

/*
def read_xyz(filename):
    """

    """

    with open(filename, "r") as f:
        frames = []
        while True:  # frame loop
            try:
                npar = int(f.readline())
                comment = f.readline()
            except ValueError:
                break
            particle_coords = []
            i = 1
            for line in f: # particle loop
                chem, x, y, z = line.split()
                particle_coords.append([x, y, z])
                i += 1
                if i > npar:
                    break
            frames.append(particle_coords)

    return asarray(frames, dtype=float)
*/

mat read_xyz(const string filename)
{
    /*
     */
    mat r(1, 3);
    return r;
}


void write_xyz(const string filename, const mat positions, const vector<string> chem_symbols, const string info, const bool append)
{
    /* 
     */
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
